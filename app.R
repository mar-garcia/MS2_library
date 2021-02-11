library(shiny)
library(Spectra)
library(CluMSID)
library(MetaboCoreUtils)
library(Rdisop)
library(MsCoreUtils)
library(CompoundDb)
library(xcms)

#db <- read.csv("database.csv")
load(list.files()[grep("MS2_library", list.files())][
  grep(".RData", list.files()[grep("MS2_library", list.files())])])

ui <- navbarPage(
  "MS2 library",
  theme = shinythemes::shinytheme("united"),  
  tabPanel("Table", 
           sidebarLayout(
             sidebarPanel(
               selectInput("polarity", label = "Polarity:", 
                           choices = list("POS" = "POS", "NEG" = "NEG"), 
                           selected = "NEG"),
               fluidRow(
                 column(3,
                        numericInput(inputId = "mz",
                                     label = "m/z value:",
                                     value = 500,
                                     step = 0.0001)),
                 
                 column(3, 
                        numericInput(inputId = "ppm",
                                     label = "ppm:",
                                     value = 1000000,
                                     step = 1)),
                 column(3, actionButton("button", "Filter by m/z"))),
               sliderInput("rt", "RT zoom:",
                           min = 0, max = 12,
                           value = c(0, 12), step = 0.1),
               sliderInput("intensity", "Relative intensity:",
                           min = 0, max = 100,
                           value = 10)
             ),
             mainPanel(
               fluidRow(DT::dataTableOutput("table")),
               fluidRow(
                 column(6, plotOutput("eic")),
                 column(6, plotOutput("ms2"))
               ),
               fluidRow(verbatimTextOutput("path"))
             )
           )
  ),
  tabPanel("Correlations",
           sidebarLayout(
             sidebarPanel(
               fluidRow(fileInput("file1", "Choose TXT file")),
               fluidRow(plotOutput("ms2_x"))),
             mainPanel(
               fluidRow(
                 column(6, DT::dataTableOutput("spectra")),
                 column(6, DT::dataTableOutput("clumsid"))),
               fluidRow(
                 column(6, plotOutput("ms2_spectra")),
                 column(6, plotOutput("ms2_clumsid"))
               )
             ))),
  navbarMenu("More",
             tabPanel("Find fragment", 
                      h3("Find spectra that contain a specific fragment ion"),
                      sidebarPanel(
                        numericInput(inputId = "frag_mz",
                                     label = "m/z value:",
                                     value = 153.0557,
                                     step = 0.0001),
                        numericInput(inputId = "frag_tol",
                                     label = "tolerance:",
                                     value = 0.001,
                                     step = 0.001)),
                      mainPanel(
                        verbatimTextOutput("frag_value")
                      )),
             tabPanel("Find loss", 
                      h3("Find spectra that contain a specific neutral loss"),
                      sidebarPanel(
                        numericInput(inputId = "loss_mz",
                                     label = "m/z value:",
                                     value = 162.0528,
                                     step = 0.0001),
                        numericInput(inputId = "loss_tol",
                                     label = "tolerance:",
                                     value = 0.001,
                                     step = 0.001)),
                      mainPanel(
                        verbatimTextOutput("loss_value")
                      ))
  ))

server <- function(input, output) {
  dbx <- reactive({
    db <- db[db$polarity == input$polarity,]
    input$button
    isolate(
      db <- db[(db$mz > input$mz - MsCoreUtils::ppm(input$mz, input$ppm)) & 
                 (db$mz < input$mz + MsCoreUtils::ppm(input$mz, input$ppm)),])
  })
  
  output$table <- DT::renderDataTable(DT::datatable({
    db <- dbx()
    db$mz <- sprintf("%.5f", round(db$mz, 5))
    db[, c("name", "formula", "adduct", "mz")]
  }, rownames= FALSE))
  
  output$eic <- renderPlot({
    db <- dbx()
    i <- input$table_rows_selected
    if (length(i) > 0){
      i <- i[length(i)]
      xdata <- Spectra(paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                       backend = MsBackendMzR())
      ms2list <- filterPrecursorMz(object = xdata, mz = db$mz[i] + 0.01 * c(-1, 1))
      rts <- rtime(ms2list)
      ms2list <- ms2list[closest(db$RT[i], rtime(ms2list))]
      
      xdata <- readMSData(files = paste0("mzML/", db$path[i], "/", db$file[i], ".mzML"), 
                          mode = "onDisk")
      chr <- chromatogram(xdata, mz = db$mz[i] + 0.01 * c(-1, 1))
      c_chr <- data.frame(
        rt = rtime(chr[[1]])/60,
        int = intensity(chr[[1]])
      )
      c_chr[is.na(c_chr)] <- 0
      plot(c_chr$rt, y = c_chr$int, type = "l", main = db$name[i], 
           xlab = "retention time", ylab = "intensity",
           xlim = c(input$rt[1], input$rt[2])) 
      points(rts/60, c_chr$int[closest(rts/60, c_chr$rt)],
             pch = 16, col = "grey")
      points(rtime(ms2list)/60, c_chr$int[closest(rtime(ms2list)/60, c_chr$rt)],
             pch = 16, col = "red")
    }
  })
  
  output$ms2 <- renderPlot({
    db <- dbx()
    i <- input$table_rows_selected
    if (length(i) > 0){
      i <- i[length(i)]
      ms2list <- spd[spd$name == db$name[i]]
      ms2list <- filterPrecursorMz(object = ms2list, mz = db$mz[i] + 0.01 * c(-1, 1))
      #ms2list <- ms2list[closest(db$RT[i], rtime(ms2list))]
      sps <- data.frame(
        mz = unlist(mz(ms2list)),
        int = unlist(intensity(ms2list))
      )
      sps$int100 <- (sps$int*100) / max(sps$int)
      idx <- which(sps$int100 >= input$intensity)
      plot(sps$mz, sps$int100, type = "h",
           xlab = "m/z", ylab = "relative intensity", 
           xlim = c(50, db$mz[i]+10),
           #xlim = c(min(sps$mz)-10, max(sps$mz)+10), 
           ylim = c(0, 110)) 
      text(sps$mz[idx], sps$int100[idx], sprintf("%.4f", round(sps$mz[idx], 4)), 
           offset = -1, pos = 2, srt = -30)
    } 
  })
  
  output$path <- renderPrint({
    db <- dbx()
    i <- input$table_rows_selected
    if (length(i) > 0){
      i <- i[length(i)]
      paste(db$path[i], db$file[i], sep = "/")
      }
  })
  
  output$spectra <- DT::renderDataTable(DT::datatable({
    req(input$file1)
    df <- read.table(input$file1$datapath)
    x_spd <- DataFrame(
      msLevel = 2L,
      polarity = 1L,
      id = "x",
      name = "x")
    x_spd$mz <- list(df[,1])
    x_spd$intensity <- list(df[,2])
    x_spd <- Spectra(x_spd)
    spdx <- filterPolarity(
      spd, 
      which(factor(c(1,2), labels = c("NEG", "POS")) == input$polarity))
    c_spd <- c(x_spd, spdx)
    tb <- cbind(c_spd$name[order(Spectra::compareSpectra(
      c_spd, ppm = 50)[,1], decreasing = T)], 
      Spectra::compareSpectra(c_spd, ppm = 50)[,1][order(
        Spectra::compareSpectra(c_spd, ppm = 50)[,1], 
        decreasing = T)]
    )[-1,]
    colnames(tb) <- c("name", "corr")
    tb[,2] <- sprintf("%.3f", round(as.numeric(tb[,2]), 3))
    return(tb)
  }))
  
  output$clumsid <- DT::renderDataTable(DT::datatable({
    req(input$file1)
    df <- read.table(input$file1$datapath)
    ms2clux <- ms2clu[lapply(ms2clu, "slot", "polarity") == input$polarity]
    
    ms2clu_i@id <- "x"
    ms2clu_i@annotation <- "x"
    ms2clu_i@precursor <- 1000
    ms2clu_i@rt <- 500
    ms2clu_i@polarity <- "x"
    ms2clu_i@spectrum <- cbind(df[,1], df[,2])
    tmp <- matrix((getSimilarities(ms2clu_i, ms2clux)[order(
      getSimilarities(ms2clu_i, ms2clux), decreasing = T)]), ncol = 1)
    tmp <- cbind(names((getSimilarities(ms2clu_i, ms2clux)[order(
      getSimilarities(ms2clu_i, ms2clux), decreasing = T)])), tmp)
    colnames(tmp) <- c("name", "corr")
    tmp[,2] <- sprintf("%.3f", round(as.numeric(tmp[,2]), 3))
    return(tmp)
  }))
  
  output$ms2_x <- renderPlot({
    req(input$file1)
    df <- read.table(input$file1$datapath)
    plot(df[,1], df[, 2], type = "h", xlab = "m/z", ylab = "intensity",
         xlim = c(min(df[,1])-10, max(df[,1]+10)),
         ylim = c(0, max(df[,2]*1.1)))
    text(df[,1], df[,2], sprintf("%.4f", round(df[,1], 4)), 
         offset = -1, pos = 2, srt = -30)
  })
  
  output$ms2_spectra <- renderPlot({
    df <- read.table(input$file1$datapath)
    x_spd <- DataFrame(
      msLevel = 2L,
      polarity = 1L,
      id = "x",
      name = "x")
    x_spd$mz <- list(df[,1])
    x_spd$intensity <- list(df[,2])
    x_spd <- Spectra(x_spd)
    spdx <- filterPolarity(
      spd, 
      which(factor(c(1,2), labels = c("NEG", "POS")) == input$polarity))
    c_spd <- c(x_spd, spdx)
    tb <- Spectra::compareSpectra(c_spd, ppm = 50)
    tb <- cbind(
      c_spd$name[order(tb[,1], decreasing = T)],
      tb[order(tb[,1], decreasing = T),1])[-1, ]
    colnames(tb) <- c("name", "corr")
    tb[,2] <- sprintf("%.3f", round(as.numeric(tb[,2]), 3))
    j <- input$spectra_rows_selected
    if (length(j) == 1){
      i <- which(spd$name == tb[j,1])
      idx <- which((unlist(intensity(spd[i]))*100)/max(unlist(intensity(spd[i]))) >= input$intensity)
      plot(unlist(mz(spd[i])), 
           (unlist(intensity(spd[i]))*100)/max(unlist(intensity(spd[i]))), 
           type = "h", main = spd$name[i],
           xlab = "m/z", ylab = "intensity", 
           xlim = c(min(unlist(mz(spd[i])))-10, max(unlist(mz(spd[i])))+10), 
           ylim = c(0, 110)
      ) 
      text(unlist(mz(spd[i]))[idx], 
           (unlist(intensity(spd[i]))[idx]*100)/max(unlist(intensity(spd[i]))), 
           sprintf("%.4f", round(unlist(mz(spd[i]))[idx], 4)), 
           offset = -1, pos = 2, srt = -30)
    }
  })
  
  output$ms2_clumsid <- renderPlot({
    df <- read.table(input$file1$datapath)
    ms2clux <- ms2clu[lapply(ms2clu, "slot", "polarity") == input$polarity]
    ms2clu_i@id <- "x"
    ms2clu_i@annotation <- "x"
    ms2clu_i@precursor <- 1000
    ms2clu_i@rt <- 500
    ms2clu_i@polarity <- "x"
    ms2clu_i@spectrum <- cbind(df[,1], df[,2])
    tmp <- getSimilarities(ms2clu_i, ms2clux)
    tmp <- cbind(names(tmp), matrix(tmp, ncol = 1))
    tmp <- tmp[order(tmp[,2], decreasing = T),]
    colnames(tmp) <- c("name", "corr")
    tmp[,2] <- sprintf("%.3f", round(as.numeric(tmp[,2]), 3))
    j <- input$clumsid_rows_selected
    names <- c()
    for(k in 1:length(ms2clux)){
      names <- c(names, ms2clux[[k]]@id)
    }
    if (length(j) == 1){
      i <- which(names == tmp[j,1])
      idx <- which(
        (ms2clux[[i]]@spectrum[,2]*100)/max(ms2clux[[i]]@spectrum[,2]) >= 
          input$intensity)
      plot(ms2clux[[i]]@spectrum[,1], 
           (ms2clux[[i]]@spectrum[,2]*100)/max(ms2clux[[i]]@spectrum[,2]), 
           type = "h", main = ms2clux[[i]]@id,
           xlab = "m/z", ylab = "intensity", 
           xlim = c(min(ms2clux[[i]]@spectrum[,1])-10, 
                    max(ms2clux[[i]]@spectrum[,1])+10), 
           ylim = c(0, 110)
      ) 
      text(ms2clux[[i]]@spectrum[idx,1], 
           (ms2clux[[i]]@spectrum[idx,2]*100)/max(ms2clux[[i]]@spectrum[,2]), 
           sprintf("%.4f", round(ms2clux[[i]]@spectrum[idx,1], 4)), 
           offset = -1, pos = 2, srt = -30)
    }
  })
  
  output$frag_value <- renderPrint({ 
    findFragment(ms2clu, mz = input$frag_mz, tolerance = input$frag_tol)
  })
  
  output$loss_value <- renderPrint({ 
    findNL(ms2clu_nl, mz = input$loss_mz, tolerance = input$loss_tol)
  })
}

shinyApp(ui = ui, server = server)