# STA2ASM Week 6: Power and sample size
# Standalone Shiny app. Requires shiny; calculations use base R.
# The displayed example R code uses MESS, but the running app does not require it.
library(shiny)

PAL <- list(teal = "#236A6E", blue = "#174F8A", rose = "#C65D7B",
            lightteal = "#EEF5F5")
fmt <- function(x, digits = 3) if (is.finite(x)) formatC(x, format = "f", digits = digits) else "Not available"
math <- function(x) withMathJax(HTML(paste0("\\(", x, "\\)")))

paired_sd <- function(sd1, sd2, rho) sqrt(sd1^2 + sd2^2 - 2 * rho * sd1 * sd2)
paired_power <- function(n, effect, sd, alpha) {
  if (n < 2 || sd <= 0 || abs(effect) < 1e-10) return(if (abs(effect) < 1e-10) alpha else NA_real_)
  tryCatch(power.t.test(n = n, delta = abs(effect), sd = sd,
                        sig.level = alpha, type = "paired", alternative = "two.sided")$power,
           error = function(e) NA_real_)
}
paired_required <- function(effect, sd, alpha, target) {
  if (abs(effect) < 1e-10) return(NA_real_)
  tryCatch(ceiling(power.t.test(power = target, delta = abs(effect), sd = sd,
                                sig.level = alpha, type = "paired", alternative = "two.sided")$n),
           error = function(e) NA_real_)
}

# MESS requires ratio >= 1, with its second group the larger group.
# The interface always defines k = n_control / n_treatment, including k < 1.
welch_arguments <- function(n1, n2, sd1, sd2) {
  reversed <- n2 < n1
  if (reversed) {
    list(n = n2, ratio = n1 / n2, sd = sd2,
         sd.ratio = sd1 / sd2, reversed = TRUE)
  } else {
    list(n = n1, ratio = n2 / n1, sd = sd1,
         sd.ratio = sd2 / sd1, reversed = FALSE)
  }
}
# Welch-Satterthwaite noncentral-t power approximation, implemented in Base R.
# The app uses this internally to run in Shinylive without MESS.
# The R code shown to students continues to use MESS::power_t_test().
welch_power <- function(n1, n2, effect, sd1, sd2, alpha) {
  if (!all(is.finite(c(n1, n2, effect, sd1, sd2, alpha))) ||
      min(n1, n2) <= 1 || min(sd1, sd2) <= 0 ||
      alpha <= 0 || alpha >= 1) return(NA_real_)
  
  v1 <- sd1^2 / n1
  v2 <- sd2^2 / n2
  se <- sqrt(v1 + v2)
  df <- (v1 + v2)^2 / (v1^2 / (n1 - 1) + v2^2 / (n2 - 1))
  critical <- qt(1 - alpha / 2, df = df)
  ncp <- abs(effect) / se
  
  # Both tails of the two-sided test (strict = TRUE in MESS).
  pt(-critical, df = df, ncp = ncp) +
    pt(critical, df = df, ncp = ncp, lower.tail = FALSE)
}

welch_sizes <- function(total, ratio) {
  n1 <- max(2L, round(total / (1 + ratio)))
  n2 <- max(2L, round(total - n1))
  c(n1 = n1, n2 = n2)
}

welch_required <- function(effect, sd1, sd2, ratio, alpha, target) {
  missing_sizes <- c(n1 = NA_real_, n2 = NA_real_)
  if (!all(is.finite(c(effect, sd1, sd2, ratio, alpha, target))) ||
      abs(effect) < 1e-10 || min(sd1, sd2, ratio) <= 0 ||
      target <= alpha || target >= 1) return(missing_sizes)
  
  # Solve for the smaller group, as in the MESS allocation convention.
  # Reverse the group SDs when the requested control/treatment ratio is < 1.
  reversed <- ratio < 1
  ratio_mess <- if (reversed) 1 / ratio else ratio
  sd_first <- if (reversed) sd2 else sd1
  sd_second <- if (reversed) sd1 else sd2
  
  power_at <- function(n) {
    welch_power(n, ratio_mess * n, effect, sd_first, sd_second, alpha)
  }
  lower <- max(2.000001, 2 / ratio_mess + 0.000001)
  upper <- max(4, lower * 2)
  while (is.finite(power_at(upper)) && power_at(upper) < target && upper < 1e7) {
    upper <- upper * 2
  }
  if (!is.finite(power_at(upper)) || power_at(upper) < target) return(missing_sizes)
  
  n_continuous <- if (power_at(lower) >= target) lower else
    tryCatch(uniroot(function(n) power_at(n) - target,
                     lower = lower, upper = upper, tol = 1e-7)$root,
             error = function(e) NA_real_)
  if (!is.finite(n_continuous)) return(missing_sizes)
  
  # Integer sample sizes: round up the smaller group and preserve the
  # allocation ratio by rounding up the larger group as well.
  first <- ceiling(n_continuous)
  second <- ceiling(first * ratio_mess)
  if (reversed) c(n1 = second, n2 = first) else c(n1 = first, n2 = second)
}

ui <- fluidPage(
  tags$head(tags$style(HTML('
    body {font-family:"Helvetica Neue",Arial,sans-serif;color:#333;background:white;}
    h1,h2,h3,h4 {color:#236A6E;font-weight:500;}
    .well {background:#EEF5F5;border:1px solid #cfdfe0;border-left:4px solid #4d8386;}
    .result {background:#EEF5F5;border-left:4px solid #236A6E;padding:12px 15px;margin:12px 0;}
    .statgrid {display:flex;flex-wrap:wrap;gap:12px;margin:10px 0;}
    .statcard {flex:1 1 145px;border:1px solid #cfdfe0;border-radius:4px;padding:10px;background:white;}
    .statcard .value {font-size:20px;font-weight:600;color:#236A6E;}
    .statcard .key {font-size:12px;color:#58696a;}
    .primary-result .statcard {background:#F4F8F8;border-left:5px solid #236A6E;}
    .primary-result .value {font-size:29px;}
    .well .control-label {font-weight:600;}
    .subtle {color:#58696a;font-size:13px;}
    .note {background:#EEF5F5;border-left:4px solid #236A6E;padding:10px 14px;margin:12px 0;}
    .btn-primary {background:#236A6E;border-color:#236A6E;}
    .shiny-input-container {max-width:100%;}
    pre {background:#fafafa;border:1px solid #dce6e6;white-space:pre-wrap;}
  '))),
  titlePanel("STA2ASM — Week 6: Power and sample size"),
  withMathJax(),
  p(class = "subtle", "Explore how study design, variability, correlation and allocation affect power and sample-size planning."),
  radioButtons("mode", "What would you like to calculate?", inline = TRUE,
               choices = c("Required sample size for target power" = "required",
                           "Power for a specified sample size" = "power"), selected = "required"),
  radioButtons("design", "Study design", inline = TRUE,
               choices = c("Paired samples" = "paired", "Independent samples (Welch)" = "welch")),
  fluidRow(
    column(3, wellPanel(
      conditionalPanel("input.design == 'paired'",
                       h4("Paired design: SBP example"),
                       numericInput("p_effect", HTML("Assumed mean paired difference \\(\\mu_D-\\mu_0\\)"), 8, step = 1),
                       numericInput("p_sd1", HTML("Baseline population SD \\(\\sigma_1\\)"), 15, min = 0.01),
                       numericInput("p_sd2", HTML("Follow-up population SD \\(\\sigma_2\\)"), 15, min = 0.01),
                       sliderInput("p_rho", HTML("Within-participant correlation \\(\\rho\\)"),
                                   min = -0.9, max = 0.9, value = 0.8, step = 0.05),
                       numericInput("p_n", HTML("Number of pairs \\(n\\)"), 30, min = 2, step = 1)
      ),
      conditionalPanel("input.design == 'welch'",
                       h4("Independent samples: hypothetical PASE study"),
                       numericInput("w_effect", HTML("Assumed mean difference \\(\\mu_1-\\mu_2\\)"), 25, step = 1),
                       numericInput("w_sd1", HTML("Treatment population SD \\(\\sigma_1\\)"), 100, min = 0.01),
                       numericInput("w_sd2", HTML("Control population SD \\(\\sigma_2\\)"), 50, min = 0.01),
                       numericInput("w_total", HTML("Total sample size \\(N=n_1+n_2\\)"), 200, min = 4, step = 1),
                       numericInput("w_ratio", HTML("Allocation ratio \\(k=n_2/n_1\\)"), 1, min = 0.1, max = 10, step = 0.1)
      ),
      selectInput("alpha", HTML("Significance level \\(\\alpha\\)"),
                  choices = c("0.01" = 0.01, "0.05" = 0.05, "0.10" = 0.10), selected = "0.05"),
      numericInput("target", HTML("Target power \\(1-\\beta\\)"), value = 0.80,
                   min = 0.1, max = 0.99, step = 0.05),
      selectInput("metric", "Graph vertical axis", choices = c("Required sample size" = "required", "Power" = "power"), selected = "power"),
      selectInput("explore", "Explore the effect of", choices = c("Sample size" = "n")),
      actionButton("reset", "Restore reading defaults")
    )),
    column(9,
           uiOutput("results"),
           h3(textOutput("plot_title")),
           plotOutput("curve", height = "385px"),
           uiOutput("explanation"),
           h4("R code for these settings (power and required sample size)"),
           tags$pre(tags$code(textOutput("rcode")))
    )
  )
)

server <- function(input, output, session) {
  # All inputs stay mounted when switching designs, so their values are preserved.
  # The graph defaults to power versus sample size, independently of
  # whether the main calculation is power or required sample size.
  observeEvent(list(input$design, input$metric), {
    choices <- if (input$design == "paired") {
      c("Number of pairs" = "n", "Mean paired difference" = "effect",
        "Baseline SD" = "sd1", "Follow-up SD" = "sd2",
        "Correlation" = "rho", "Significance level" = "alpha")
    } else {
      c("Total sample size" = "n", "Mean difference" = "effect",
        "Treatment SD" = "sd1", "Control SD" = "sd2",
        "Allocation ratio" = "ratio", "Significance level" = "alpha")
    }
    if (input$metric == "required") {
      choices <- c(choices[names(choices) != "Number of pairs" &
                             names(choices) != "Total sample size"],
                   "Target power" = "target")
    }
    updateSelectInput(session, "explore", choices = choices, selected = if (input$metric == "power") "n" else "effect")
  }, ignoreInit = FALSE)
  observeEvent(input$reset, {
    updateRadioButtons(session, "mode", selected = "required")
    updateSelectInput(session, "metric", selected = "power")
    updateSelectInput(session, "alpha", selected = "0.05")
    updateNumericInput(session, "target", value = 0.8)
    updateSelectInput(session, "explore", selected = "n")
    if (input$design == "paired") {
      updateNumericInput(session, "p_effect", value = 8)
      updateNumericInput(session, "p_sd1", value = 15)
      updateNumericInput(session, "p_sd2", value = 15)
      updateSliderInput(session, "p_rho", value = 0.8)
      updateNumericInput(session, "p_n", value = 30)
    } else {
      updateNumericInput(session, "w_effect", value = 25)
      updateNumericInput(session, "w_sd1", value = 100)
      updateNumericInput(session, "w_sd2", value = 50)
      updateNumericInput(session, "w_total", value = 200)
      updateNumericInput(session, "w_ratio", value = 1)
    }
  })
  settings <- reactive({
    alpha <- as.numeric(input$alpha)
    req(is.finite(alpha), alpha > 0, alpha < 1,
        is.finite(input$target), input$target > alpha, input$target < 1)
    if (input$design == "paired") {
      req(input$p_sd1 > 0, input$p_sd2 > 0, input$p_n >= 2,
          abs(input$p_rho) < 1, is.finite(input$p_effect))
      list(design = "paired", effect = input$p_effect, sd1 = input$p_sd1,
           sd2 = input$p_sd2, rho = input$p_rho, n = as.integer(input$p_n),
           alpha = alpha, target = input$target)
    } else {
      req(input$w_sd1 > 0, input$w_sd2 > 0, input$w_total >= 4,
          input$w_ratio > 0, is.finite(input$w_effect))
      list(design = "welch", effect = input$w_effect, sd1 = input$w_sd1,
           sd2 = input$w_sd2, n = as.integer(input$w_total),
           ratio = input$w_ratio, alpha = alpha, target = input$target)
    }
  })
  calculation <- reactive({
    d <- settings()
    if (d$design == "paired") {
      sdD <- paired_sd(d$sd1, d$sd2, d$rho)
      list(power = paired_power(d$n, d$effect, sdD, d$alpha),
           required = paired_required(d$effect, sdD, d$alpha, d$target),
           sdD = sdD, se = sdD / sqrt(d$n))
    } else {
      ns <- welch_sizes(d$n, d$ratio)
      list(power = welch_power(ns[1], ns[2], d$effect, d$sd1, d$sd2, d$alpha),
           required = welch_required(d$effect, d$sd1, d$sd2, d$ratio, d$alpha, d$target),
           sizes = ns, se = sqrt(d$sd1^2 / ns[1] + d$sd2^2 / ns[2]),
           optimum = d$sd2 / d$sd1)
    }
  })
  card <- function(label, value) div(class = "statcard", div(class = "key", label),
                                     div(class = "value", value))
  output$results <- renderUI({
    d <- settings(); r <- calculation()
    if (d$design == "paired") {
      main <- if (input$mode == "required")
        card("REQUIRED NUMBER OF PAIRS", if (is.finite(r$required)) r$required else "Not available")
      else card("CALCULATED POWER", fmt(r$power))
      div(class = "result",
          div(class = "primary-result", main),
          div(class = "statgrid",
              card("SD of paired differences", fmt(r$sdD, 2)),
              card("Standard error at selected n", fmt(r$se, 2)),
              if (input$mode == "required") card("Target power", fmt(d$target, 2))
              else card("Number of pairs", d$n)))
    } else {
      main <- if (input$mode == "required") {
        if (all(is.finite(r$required))) {
          div(class = "statgrid",
              card("REQUIRED TREATMENT PARTICIPANTS", r$required[1]),
              card("REQUIRED CONTROL PARTICIPANTS", r$required[2]),
              card("REQUIRED TOTAL", sum(r$required)))
        } else card("REQUIRED SAMPLE SIZES", "Not available")
      } else card("CALCULATED POWER", fmt(r$power))
      div(class = "result",
          div(class = "primary-result", main),
          div(class = "statgrid",
              card("Current group sizes (treatment / control)", paste(r$sizes[1], "/", r$sizes[2])),
              card("Standard error at selected sizes", fmt(r$se, 2)),
              card("Optimal allocation ratio", fmt(r$optimum, 2)),
              if (input$mode == "required") card("Target power", fmt(d$target, 2))))
    }
  })
  output$plot_title <- renderText({
    axis_name <- switch(input$explore,
                        n = "sample size", effect = "assumed mean difference",
                        sd1 = "population SD 1", sd2 = "population SD 2",
                        rho = "correlation", ratio = "allocation ratio",
                        alpha = "significance level", target = "target power", "parameter")
    paste(if (input$metric == "required") "Required sample size" else "Power",
          "versus", axis_name)
  })
  curve_data <- reactive({
    d <- settings(); variable <- input$explore
    req(variable)
    current <- switch(variable, n = d$n, effect = abs(d$effect), sd1 = d$sd1,
                      sd2 = d$sd2, rho = d$rho, ratio = d$ratio, alpha = d$alpha, target = d$target)
    required_n <- if (d$design == "paired") calculation()$required else sum(calculation()$required)
    n_max <- max(c(80, d$n * 2, if (is.finite(required_n)) required_n * 1.35), na.rm = TRUE)
    xx <- switch(variable,
                 n = seq(if (d$design == "paired") 2 else 8, n_max, length.out = 75),
                 effect = seq(0, max(2, abs(d$effect) * 2), length.out = 75),
                 sd1 = seq(max(0.1, d$sd1 * 0.2), d$sd1 * 2.2, length.out = 75),
                 sd2 = seq(max(0.1, d$sd2 * 0.2), d$sd2 * 2.2, length.out = 75),
                 rho = seq(-0.9, 0.9, length.out = 75),
                 ratio = seq(0.2, max(3, d$ratio * 1.5), length.out = 75),
                 alpha = seq(0.005, 0.2, length.out = 75),
                 target = seq(max(d$alpha + 0.02, 0.15), 0.95, length.out = 75))
    yy <- vapply(xx, function(z) {
      a <- d
      if (variable == "n") a$n <- round(z)
      if (variable == "effect") a$effect <- z
      if (variable == "sd1") a$sd1 <- z
      if (variable == "sd2") a$sd2 <- z
      if (variable == "rho") a$rho <- z
      if (variable == "ratio") a$ratio <- z
      if (variable == "alpha") a$alpha <- z
      if (variable == "target") a$target <- z
      if (input$metric == "required") {
        if (a$target <= a$alpha) return(NA_real_)
        if (a$design == "paired") {
          paired_required(a$effect, paired_sd(a$sd1, a$sd2, a$rho), a$alpha, a$target)
        } else {
          sum(welch_required(a$effect, a$sd1, a$sd2, a$ratio, a$alpha, a$target))
        }
      } else if (a$design == "paired") {
        paired_power(a$n, a$effect, paired_sd(a$sd1, a$sd2, a$rho), a$alpha)
      } else {
        ns <- welch_sizes(a$n, a$ratio)
        welch_power(ns[1], ns[2], a$effect, a$sd1, a$sd2, a$alpha)
      }
    }, numeric(1))
    list(x = xx, y = yy, current = current)
  })
  output$curve <- renderPlot({
    d <- settings(); r <- calculation(); curve <- curve_data()
    old <- par(mar = c(4.5, 4.7, 1.1, 1.5)); on.exit(par(old))
    labels <- c(n = if (d$design == "paired") "Number of pairs" else "Total sample size",
                effect = "Absolute assumed mean difference", sd1 = "Population SD 1",
                sd2 = "Population SD 2", rho = "Within-participant correlation",
                ratio = "Allocation ratio (control / treatment)", alpha = "Significance level",
                target = "Target power")
    is_required <- input$metric == "required"
    y_top <- if (is_required) max(5, curve$y[is.finite(curve$y)], na.rm = TRUE) * 1.1 else 1
    plot(curve$x, curve$y, type = "l", col = PAL$teal, lwd = 3,
         xlab = labels[[input$explore]],
         ylab = if (is_required) if (d$design == "paired") "Required number of pairs" else "Required total sample size" else "Power",
         ylim = c(0, y_top),
         las = 1, bty = "l")
    if (!is_required) abline(h = d$target, col = PAL$blue, lty = 2, lwd = 2)
    # In sample-size mode, mark the required sample size on a power-vs-size plot.
    mark_required <- input$mode == "required" && !is_required && input$explore == "n"
    marker_x <- if (mark_required) {
      if (d$design == "paired") r$required else sum(r$required)
    } else curve$current
    if (is.finite(marker_x)) abline(v = marker_x, col = PAL$rose, lty = 3, lwd = 2)
    selected_y <- if (is_required) {
      if (d$design == "paired") r$required else sum(r$required)
    } else r$power
    if (mark_required) {
      selected_y <- if (d$design == "paired") {
        paired_power(r$required, d$effect, r$sdD, d$alpha)
      } else {
        welch_power(r$required[1], r$required[2], d$effect, d$sd1, d$sd2, d$alpha)
      }
    }
    if (is.finite(marker_x) && is.finite(selected_y))
      points(marker_x, selected_y, pch = 21, bg = PAL$teal, cex = 1.5)
    if (d$design == "welch" && input$explore == "ratio" &&
        r$optimum >= min(curve$x) && r$optimum <= max(curve$x))
      abline(v = r$optimum, col = PAL$blue, lty = 2, lwd = 2)
    legend_labels <- c(if (is_required) "Required sample size" else "Power",
                       if (!is_required) "Target power", if (mark_required) "Required sample size" else "Current setting",
                       if (d$design == "welch" && input$explore == "ratio") "Optimal allocation")
    legend_colors <- c(PAL$teal, if (!is_required) PAL$blue, PAL$rose,
                       if (d$design == "welch" && input$explore == "ratio") PAL$blue)
    legend("bottomright", inset = c(0.02, 0.02),
           legend = legend_labels, col = legend_colors,
           lty = c(1, if (!is_required) 2, 3,
                   if (d$design == "welch" && input$explore == "ratio") 2),
           lwd = 2, bty = "n", cex = 0.83)
  })
  
  output$explanation <- renderUI({
    d <- settings()
    if (d$design == "paired") {
      div(class = "note", withMathJax(HTML(
        "Holding \\(\\sigma_1\\) and \\(\\sigma_2\\) fixed, increasing the correlation \\(\\rho\\) reduces \\(\\sigma_D\\), which increases power. For example, try comparing \\(\\rho = 0.8\\) with \\(\\rho = -0.8\\)."
      )))
    } else {
      div(class = "note", withMathJax(HTML(
        "When the population standard deviations differ between groups, allocating more participants to the group with greater variability can increase power, even when the total sample size remains unchanged. With equal recruitment costs, the allocation ratio that minimises the standard error is \\(k = \\sigma_2 / \\sigma_1\\)."
      )))
    }
  })
  
  output$rcode <- renderText({
    d <- settings()
    if (d$design == "paired") {
      paste0("# Paired power: assumed population values\n",
             "sd1 <- ", d$sd1, "\nsd2 <- ", d$sd2, "\nrho <- ", d$rho,
             "\nsigma_D <- sqrt(sd1^2 + sd2^2 - 2 * rho * sd1 * sd2)\n",
             "power.t.test(n = ", d$n, ", delta = ", abs(d$effect),
             ", sd = sigma_D, sig.level = ", d$alpha,
             ", type = \"paired\", alternative = \"two.sided\")\n",
             "power.t.test(power = ", d$target, ", delta = ", abs(d$effect),
             ", sd = sigma_D, sig.level = ", d$alpha,
             ", type = \"paired\", alternative = \"two.sided\")")
    } else {
      r <- calculation()
      a <- welch_arguments(r$sizes[1], r$sizes[2], d$sd1, d$sd2)
      note <- if (a$reversed) {
        paste0("# Group order reversed automatically because the selected k = ",
               d$ratio, " is below 1. MESS requires ratio >= 1.\n",
               "# MESS group 1 = Control; MESS group 2 = Treatment.\n")
      } else {
        "# MESS group 1 = Treatment; MESS group 2 = Control.\n"
      }
      req_a <- if (d$ratio < 1) list(sd = d$sd2, sd.ratio = d$sd1 / d$sd2, ratio = 1 / d$ratio) else
        list(sd = d$sd1, sd.ratio = d$sd2 / d$sd1, ratio = d$ratio)
      required_note <- if (d$ratio < 1) "# For sample-size planning, reverse the groups because the chosen allocation ratio is below 1.\n" else ""
      paste0("# Run this code in R with the MESS package installed.\n",
             "library(MESS)\n",
             "# Selected group sizes: Treatment = ", r$sizes[1],
             ", Control = ", r$sizes[2], "\n", note,
             "# Power for the current allocation\n",
             "power_t_test(n = ", a$n, ", delta = ", abs(d$effect),
             ", sd = ", a$sd, ", sd.ratio = ", signif(a$sd.ratio, 8),
             ",\n             ratio = ", signif(a$ratio, 8),
             ", sig.level = ", d$alpha, ", power = NULL,\n",
             "             type = \"two.sample\", alternative = \"two.sided\",\n",
             "             df.method = \"welch\", strict = TRUE)\n\n",
             "# Required sample sizes for target power ", d$target, "\n", required_note,
             "power_t_test(n = NULL, delta = ", abs(d$effect),
             ", sd = ", req_a$sd, ", sd.ratio = ", signif(req_a$sd.ratio, 8),
             ",\n             ratio = ", signif(req_a$ratio, 8),
             ", sig.level = ", d$alpha, ", power = ", d$target, ",\n",
             "             type = \"two.sample\", alternative = \"two.sided\",\n",
             "             df.method = \"welch\", strict = TRUE)\n",
             "# Round required sample sizes up to whole participants.")
    }
  })
}

shinyApp(ui, server)
