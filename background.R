library(shiny)

#The path is to the folder, not the app file (app should be called app.R or server.R/ui.R)
runApp("app.R", launch.browser = T, host = "127.0.0.1", port = 8100)
