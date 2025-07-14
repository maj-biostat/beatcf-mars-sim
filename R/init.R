# load libs, load config, initialise log

source("R/libs.R")

is_html <- knitr::is_html_output()
# is_pdf <- knitr::is_latex_output()
# is_word <- !is_html & !is_pdf

ggplot2::theme_set(ggplot2::theme_bw())
ggplot2::theme_update(text = element_text(size = 8))
ggplot2::theme_update(legend.position = "bottom")
# ggplot2::theme_update(legend.title = element_blank())
ggplot2::theme_update(axis.text.x = element_text(size = 8))
ggplot2::theme_update(axis.text.y = element_text(size = 8))

# g_cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Logs
f_log <- file.path("./logs", "log.txt")
log_appender(appender_file(f_log))
# message(Sys.time(), " Log file initialised ", f_log)
log_info("*** START UP ***")

