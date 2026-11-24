# Australian Wine Sales Explorer
# Keep AustralianWines.csv in the same folder as this file.
library(shiny)
library(fpp3)
library(DT)

wine_raw <- readr::read_csv(
  "AustralianWines.csv", na = c("", "NA", "*"),
  show_col_types = FALSE
)
names(wine_raw) <- trimws(names(wine_raw))
# Parse the source's English month abbreviations without a two-digit year cutoff.
wine <- wine_raw |>
  pivot_longer(-Month, names_to = "Varietal", values_to = "Sales") |>
  mutate(Month = yearmonth(as.Date(sprintf(
    "%04d-%02d-01", 1900L + as.integer(substr(Month, 5, 6)),
    match(substr(Month, 1, 3), month.abb)
  )))) |>
  as_tsibble(index = Month, key = Varietal)
stopifnot(!anyNA(wine$Month), is.numeric(wine$Sales))

# Impute only the window being fitted. Validation actuals are never imputed.
prepare_fit_data <- function(data) {
  data |>
    group_by(Varietal) |>
    mutate(Sales = {
      observed <- which(is.finite(Sales))
      if (length(observed) < 2) stop("Not enough observed sales to fit models.")
      approx(observed, Sales[observed], seq_along(Sales), rule = 2)$y
    }) |>
    ungroup()
}

fit_models <- function(data) {
  prepare_fit_data(data) |>
    model(
      TSLM = TSLM(Sales ~ trend() + season()),
      ETS = ETS(Sales),
      ARIMA = ARIMA(Sales)
    )
}

# Use original, observed sales for scoring, excluding imputed targets.
score_predictions <- function(predictions, actual, window) {
  predictions |>
    as_tibble() |>
    select(Varietal, Month, .model, Prediction) |>
    left_join(as_tibble(actual), by = c("Varietal", "Month")) |>
    group_by(Varietal, .model) |>
    summarise(
      N = sum(is.finite(Sales) & is.finite(Prediction)),
      RMSE = if (N > 0) sqrt(mean((Sales - Prediction)^2, na.rm = TRUE)) else NA_real_,
      MAE = if (N > 0) mean(abs(Sales - Prediction), na.rm = TRUE) else NA_real_,
      MAPE = if (any(is.finite(Sales) & Sales != 0 & is.finite(Prediction)))
        mean(abs((Sales - Prediction) / Sales)[is.finite(Sales) & Sales != 0 & is.finite(Prediction)]) * 100
        else NA_real_,
      .groups = "drop"
    ) |>
    rename(Model = .model) |>
    mutate(Window = window, .after = Model) |>
    arrange(Varietal, RMSE)
}

specifications <- function(fits) {
  as_tibble(fits) |>
    mutate(across(-Varietal, ~ purrr::map_chr(.x, function(m) format(m)))) |>
    pivot_longer(-Varietal, names_to = "Model", values_to = "Specification") |>
    mutate(Specification = if_else(Model == "TSLM",
      "Linear trend + monthly seasonal effects", Specification))
}

forecast_chart <- function(fc, actual, boundary, title) {
  bands <- as_tibble(fc) |>
    mutate(Date = as.Date(Month),
           Lower = as.numeric(quantile(Sales, 0.025)),
           Upper = as.numeric(quantile(Sales, 0.975)))
  actual <- actual |>
    as_tibble() |>
    mutate(Date = as.Date(Month))
  ggplot() +
    geom_ribbon(data = bands,
      aes(Date, ymin = Lower, ymax = Upper, fill = "95% prediction interval"), alpha = 0.25) +
    geom_line(data = actual, aes(Date, Sales, colour = "Observed sales"), na.rm = TRUE) +
    geom_line(data = bands, aes(Date, .mean, colour = "Forecast mean"), linewidth = 0.7) +
    geom_vline(xintercept = as.Date(boundary), linetype = "dashed", colour = "grey45") +
    facet_grid(Varietal ~ .model, scales = "free_y") +
    scale_colour_manual(values = c("Observed sales" = "#303030", "Forecast mean" = "#0072B2")) +
    scale_fill_manual(values = c("95% prediction interval" = "#56B4E9")) +
    labs(title = title, subtitle = "Dashed line marks the first forecast month.",
         x = NULL, y = "Sales (source units)", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
}

metric_table <- function(data) {
  datatable(data, rownames = FALSE,
    options = list(dom = "t", pageLength = 36, scrollX = TRUE)) |>
    formatRound(c("RMSE", "MAE", "MAPE"), 2)
}

ui <- fluidPage(
  titlePanel("Australian Wine Sales Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("varietals", "Wine varieties:", sort(unique(wine$Varietal)),
                  selected = c("Red", "sparkling"), multiple = TRUE),
      dateRangeInput("dates", "Data window (whole months):",
        start = as.Date(min(wine$Month)), end = as.Date(max(wine$Month)),
        min = as.Date(min(wine$Month)), max = as.Date(max(wine$Month))),
      checkboxGroupInput("models", "Models to compare:",
        choices = c("TSLM", "ETS", "ARIMA"), selected = c("TSLM", "ETS", "ARIMA")),
      sliderInput("validation_months", "Months held out for validation:", 6, 36, 12, step = 1),
      sliderInput("horizon", "Future forecast horizon (months):", 1, 36, 12, step = 1),
      helpText("Keep at least 36 training months. Validation uses the last months of the data window."),
      helpText("Future forecasts refit models on the entire selected window and start after its last month."),
      downloadButton("download_metrics", "Download accuracy CSV")
    ),
    mainPanel(
      verbatimTextOutput("split_summary"),
      textOutput("missing_summary"),
      tabsetPanel(
        tabPanel("Overview", plotOutput("sales_plot", height = "450px")),
        tabPanel("Compare models",
          h3("Best validation RMSE among selected models"),
          DTOutput("best_models"),
          helpText("A winner applies to this validation window; it is not a guarantee of future performance. Exact ties are retained."),
          h3("Training accuracy"), DTOutput("training_accuracy"),
          h3("Validation accuracy"), DTOutput("validation_accuracy"),
          helpText("Lower is better. N counts scored observations; MAPE is a percentage and excludes zero actuals. Missing actuals are not scored."),
          h3("Specifications fitted on training data"), DTOutput("model_specs"),
          p("ETS components are error, trend, seasonality: A = additive, M = multiplicative, N = none; Ad = damped additive trend."),
          p("ARIMA(p,d,q)(P,D,Q)[12] shows nonseasonal and seasonal orders; [12] denotes monthly seasonality. TSLM uses a linear trend and monthly effects.")
        ),
        tabPanel("Validation forecasts", uiOutput("validation_plot_ui")),
        tabPanel("Future forecasts",
          p("These models are refitted on all selected data, including the former validation months. Their specifications can differ from the training fits."),
          uiOutput("future_plot_ui"),
          h3("Specifications refitted on the full window"), DTOutput("future_specs")
        ),
        tabPanel("About",
          h3("Data and method"),
          p("Source: AustralianWines.csv supplied for ADS 506. Monthly observations cover January 1980 to December 1994. Sales units are retained as supplied."),
          p("The source uses * for missing Rose sales in July and August 1994. For fitting only, missing values are linearly interpolated within each variety and fitting window; endpoint gaps use the nearest observed value. Validation targets are never filled."),
          p("ETS and ARIMA use automatic specification selection. TSLM uses a specified trend-plus-season model. Training metrics use fitted values; validation forecasts are made once from the training endpoint, without updating from held-out observations."),
          p("Training metrics exclude imputed targets. Forecast intervals are model-based 95% prediction intervals and do not account for imputation or model-selection uncertainty. Negative lower limits are retained rather than silently clipped."),
          p("Date controls include whole months. Future forecasts start after the selected data window, not necessarily after December 1994."),
          h3("Reproduce the initial comparison"),
          p("Select Red and sparkling, January 1980 through December 1994, all three models, 12 validation months, and a 12-month future horizon."),
          h3("Acknowledgments"),
          p("App code drafted and revised with OpenAI ChatGPT/Codex assistance. The author is responsible for reviewing the methods and interpreting the results."),
          tags$a(href = "https://fable.tidyverts.org/", "fable documentation"), tags$br(),
          tags$a(href = "https://otexts.com/fpp3/", "Forecasting: Principles and Practice (3rd edition)")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  filtered_wine <- reactive({
    req(input$dates)
    validate(need(length(input$varietals) > 0, "Select at least one wine variety."),
             need(!anyNA(input$dates), "Choose both dates."),
             need(input$dates[1] <= input$dates[2], "The start date must precede the end date."))
    data <- wine |> filter(Varietal %in% input$varietals,
      Month >= yearmonth(input$dates[1]), Month <= yearmonth(input$dates[2]))
    validate(need(nrow(data) > 0, "Choose a window containing observations."))
    data
  })
  split_data <- reactive({
    data <- filtered_wine()
    months <- sort(unique(data$Month))
    validate(need(length(months) >= input$validation_months + 36,
      "Choose a longer date range or fewer validation months: at least 36 training months are required."))
    cutoff <- months[length(months) - input$validation_months]
    list(train = filter(data, Month <= cutoff), validation = filter(data, Month > cutoff))
  })
  selected_models <- reactive({
    validate(need(length(input$models) > 0, "Select at least one forecasting model."))
    input$models
  })
  # Fit once per data window; changing displayed models reuses these results.
  all_training_models <- reactive({
    withProgress(message = "Fitting training models...", value = 0.5,
                 fit_models(split_data()$train))
  })
  training_models <- reactive({
    chosen <- selected_models()
    select(all_training_models(), Varietal, all_of(chosen))
  })
  all_future_models <- reactive({
    split_data() # enforce the same minimum-history rule for all tabs
    withProgress(message = "Refitting on the full window...", value = 0.5,
                 fit_models(filtered_wine()))
  })
  future_models <- reactive({
    chosen <- selected_models()
    select(all_future_models(), Varietal, all_of(chosen))
  })
  validation_forecasts <- reactive({
    forecast(training_models(), h = input$validation_months)
  })
  future_forecasts <- reactive({forecast(future_models(), h = input$horizon)})
  training_metrics <- reactive({
    predictions <- augment(training_models()) |> rename(Prediction = .fitted)
    score_predictions(predictions, split_data()$train, "Training")
  })
  validation_metrics <- reactive({
    predictions <- as_tibble(validation_forecasts()) |> rename(Prediction = .mean)
    score_predictions(predictions, split_data()$validation, "Validation")
  })
  output$split_summary <- renderText({
    parts <- split_data()
    paste0("Training: ", min(parts$train$Month), " through ", max(parts$train$Month),
      "\nValidation: ", min(parts$validation$Month), " through ", max(parts$validation$Month),
      "\nFuture forecast: ", max(filtered_wine()$Month) + 1, " through ",
      max(filtered_wine()$Month) + input$horizon)
  })
  output$missing_summary <- renderText({
    paste("Missing sales in selected window:", sum(is.na(filtered_wine()$Sales)),
          "- shown as gaps and excluded from accuracy targets.")
  })
  output$sales_plot <- renderPlot({
    ggplot(filtered_wine(), aes(as.Date(Month), Sales, colour = Varietal)) +
      geom_line(na.rm = TRUE) +
      labs(title = "Monthly wine sales", x = NULL, y = "Sales (source units)", colour = "Wine variety") +
      theme_minimal(base_size = 13) + theme(legend.position = "bottom")
  })
  output$training_accuracy <- renderDT({metric_table(training_metrics())})
  output$validation_accuracy <- renderDT({metric_table(validation_metrics())})
  output$best_models <- renderDT({
    winners <- validation_metrics() |> filter(is.finite(RMSE), N > 0) |>
      group_by(Varietal) |> slice_min(RMSE, n = 1, with_ties = TRUE) |> ungroup()
    validate(need(nrow(winners) > 0, "No valid scores are available for these settings."))
    metric_table(winners)
  })
  output$model_specs <- renderDT({
    datatable(specifications(training_models()), rownames = FALSE,
              options = list(dom = "t", pageLength = 18, scrollX = TRUE))
  })
  output$future_specs <- renderDT({
    datatable(specifications(future_models()), rownames = FALSE,
              options = list(dom = "t", pageLength = 18, scrollX = TRUE))
  })
  output$validation_plot_ui <- renderUI({
    plotOutput("validation_plot", height = paste0(max(400, 260 * length(input$varietals)), "px"))
  })
  output$future_plot_ui <- renderUI({
    plotOutput("future_plot", height = paste0(max(400, 260 * length(input$varietals)), "px"))
  })
  output$validation_plot <- renderPlot({
    boundary <- min(split_data()$validation$Month)
    recent <- filtered_wine() |> filter(Month >= boundary - 24)
    forecast_chart(validation_forecasts(), recent, boundary, "Validation forecasts versus actual sales")
  })
  output$future_plot <- renderPlot({
    boundary <- max(filtered_wine()$Month) + 1
    recent <- filtered_wine() |> filter(Month >= boundary - 36)
    forecast_chart(future_forecasts(), recent, boundary, "Future forecasts from the full selected window")
  })
  output$download_metrics <- downloadHandler(
    filename = function() "wine-model-accuracy.csv",
    content = function(file) {
      parts <- split_data()
      data <- bind_rows(training_metrics(), validation_metrics()) |>
        mutate(Training_start = as.character(min(parts$train$Month)),
          Training_end = as.character(max(parts$train$Month)),
          Validation_start = as.character(min(parts$validation$Month)),
          Validation_end = as.character(max(parts$validation$Month)))
      readr::write_csv(data, file)
    }
  )
}

shinyApp(ui, server)
