library(shiny)
library(Spectra)
library(CluMSID)
library(MetaboCoreUtils)
library(Rdisop)
library(MsCoreUtils)
library(CompoundDb)
library(xcms)

#db <- read.csv("database.csv")
load("MS2_library.RData")
rm(startpoint)

ui <- navbarPage("MS2 library",
                 theme = shinythemes::shinytheme("united"),  
                 tabPanel("Table", 
                          fluidRow(DT::dataTableOutput("table")),
                          fluidRow(
                            column(6, plotOutput("eic")),
                            column(6, plotOutput("ms2"))
                          )),
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
                            tabPanel("Summary", "Summary tab contents..."),
                            tabPanel("Table", "Table tab contents...")
                 ))

server <- function(input, output) {
  output$table <- DT::renderDataTable(DT::datatable({
    db[, c("name", "formula")]
  }, rownames= FALSE))
  
  output$eic <- renderPlot({
    i <- input$table_rows_selected
    if (length(i) == 1){
      c_mz <- unlist(mass2mz(getMolecule(db$formula[i])$exactmass, db$adduct[i]))
      
      xdata <- Spectra(paste0("mzML/", db$file[i], ".mzML"), 
                       backend = MsBackendMzR())
      ms2list <- filterPrecursorMz(object = xdata, mz = c_mz + 0.01 * c(-1, 1))
      rts <- rtime(ms2list)
      ms2list <- ms2list[closest(db$RT[i], rtime(ms2list))]
      
      xdata <- readMSData(files = paste0("mzML/", db$file[i], ".mzML"), 
                          mode = "onDisk")
      chr <- chromatogram(xdata, mz = c_mz + 0.01 * c(-1, 1))
      c_chr <- data.frame(
        rt = rtime(chr[[1]])/60,
        int = intensity(chr[[1]])
      )
      c_chr[is.na(c_chr)] <- 0
      plot(c_chr$rt, y = c_chr$int, type = "l", 
           xlab = "retention time", ylab = "intensity", main = db$name[i]) 
      points(rts/60, c_chr$int[closest(rts/60, c_chr$rt)],
             pch = 16, col = "grey")
      points(rtime(ms2list)/60, c_chr$int[closest(rtime(ms2list)/60, c_chr$rt)],
             pch = 16, col = "red")
    }
  })
  
  output$ms2 <- renderPlot({
    i <- input$table_rows_selected
    if (length(i) == 1){
      c_mz <- unlist(mass2mz(getMolecule(db$formula[i])$exactmass, db$adduct[i]))
      xdata <- Spectra(paste0("mzML/", db$file[i], ".mzML"), 
                       backend = MsBackendMzR())
      ms2list <- filterPrecursorMz(object = xdata, mz = c_mz + 0.01 * c(-1, 1))
      rts <- rtime(ms2list)
      ms2list <- ms2list[closest(db$RT[i], rtime(ms2list))]
      sps <- data.frame(
        mz = unlist(mz(ms2list)),
        int = unlist(intensity(ms2list))
      )
      sps$int100 <- (sps$int*100) / max(sps$int)
      c_frag <- c(db$formula[i], unlist(strsplit(db$fragments[i], "; ")))
      c_frag_mz <- c()
      for(j in 1:length(c_frag)){
        c_frag_mz <- c(c_frag_mz, unlist(mass2mz(getMolecule(c_frag[j])$exactmass, "[M-H]-")))
      }
      idx <- unlist(matchWithPpm(c_frag_mz, sps$mz, ppm = 10))
      plot(sps$mz[idx], sps$int100[idx], type = "h",
           xlab = "m/z", ylab = "relative intensity", 
           xlim = c(min(sps$mz), max(sps$mz)), ylim = c(0, 110)) 
      text(sps$mz[idx], sps$int100[idx], sprintf("%.4f", round(sps$mz[idx], 4)), 
           offset = -1, pos = 2, srt = -30)
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
    c_spd <- c(x_spd, spd)
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
    ms2clu_i@id <- "x"
    ms2clu_i@annotation <- "x"
    ms2clu_i@precursor <- 1000
    ms2clu_i@rt <- 500
    ms2clu_i@polarity <- "x"
    ms2clu_i@spectrum <- cbind(df[,1], df[,2])
    tmp <- matrix((getSimilarities(ms2clu_i, ms2clu)[order(
      getSimilarities(ms2clu_i, ms2clu), decreasing = T)]), ncol = 1)
    tmp <- cbind(names((getSimilarities(ms2clu_i, ms2clu)[order(
      getSimilarities(ms2clu_i, ms2clu), decreasing = T)])), tmp)
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
    j <- input$spectra_rows_selected
    if (length(j) == 1){
      i <- which(spd$name == tb[j,1])
      plot(unlist(mz(spd[i])), unlist(intensity(spd[i])), 
           type = "h", main = db$name[i],
           xlab = "m/z", ylab = "intensity", 
           xlim = c(min(unlist(mz(spd[i])))-10, max(unlist(mz(spd[i])))+10), 
           ylim = c(0, max(unlist(intensity(spd[i])*1.1)))
           ) 
      text(unlist(mz(spd[i])), unlist(intensity(spd[i])), 
           sprintf("%.4f", round(unlist(mz(spd[i])), 4)), 
           offset = -1, pos = 2, srt = -30)
    }
    })
  
  output$ms2_clumsid <- renderPlot({
    i <- input$clumsid_rows_selected
    names <- c()
    for(k in 1:length(ms2clu)){
      names <- c(names, ms2clu[[k]]@id)
    }
    if (length(i) == 1){
      i <- which(names == tb[j,1])
      plot(ms2clu[[i]]@spectrum[,1], ms2clu[[i]]@spectrum[,2], 
           type = "h", main = ms2clu[[i]]@id,
           xlab = "m/z", ylab = "intensity", 
           xlim = c(min(ms2clu[[i]]@spectrum[,1])-10, max(ms2clu[[i]]@spectrum[,1])+10), 
           ylim = c(0, max(ms2clu[[i]]@spectrum[,2]*1.1))
      ) 
      text(ms2clu[[i]]@spectrum[,1], ms2clu[[i]]@spectrum[,2], 
           sprintf("%.4f", round(ms2clu[[i]]@spectrum[,1], 4)), 
           offset = -1, pos = 2, srt = -30)
    }
  })
}

shinyApp(ui = ui, server = server)