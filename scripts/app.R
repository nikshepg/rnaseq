library(shiny)
library(shinythemes)
library(colourpicker)
library(scales)

library(rsconnect)

plotOntology <- function(a, low_color, high_color, font) {
  ggplot(a[1:10, ], aes(y = fct_inorder(Term), x = Combined.Score, fill = Combined.Score)) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient(low= low_color,high = high_color) +
    scale_y_discrete(labels = label_wrap(40)) +
    labs(x = "Combined Score", y = "Process") +
    theme(text = element_text(size = 12, family = font))
}

MDSplot <- function(i) {

  if (i == 1) {
    lcpm <- cpm(x, log=TRUE)
    par(mfrow=c(1,2))
    col.group <-group

    levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
    col.group <- as.character(col.group)
    plotMDS(lcpm, labels=c("T1S1","T1S2","T1S3","T1S4","T1S5","T1S6","T1S7","T1S8","T2S1","T2S2","T2S3","T2S4","T2S5","T2S6"), col=col.group)
    title(main="Sample groups")
  }

  else if (i == 2) {
    boxplot(lcpm, las=2, col=col, main="")
    title(main="Normalised Data", ylab="Log-cpm")
  }

  else if (i == 3) {
    lcpm.cutoff <- log2(10/M + 2/L)
    nsamples <- ncol(x1)
    col <- brewer.pal(nsamples, "Paired")
    par(mfrow=c(1,2))
    plot(density(lcpm1[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(main="Raw data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
      den <- density(lcpm1[,i])
      lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", samplenames, text.col=col, bty="n")
  }
  else if(i ==4){
    lcpm <- cpm(x, log=TRUE)
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(main="Filtered data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
      den <- density(lcpm[,i])
      lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", samplenames, text.col=col, bty="n")
  }}

DEGgraph <- function(i, p) {
  dt <- decideTests(tfit, p.value = p)
  n.value <- length(which(FirstVsSecond1$adj.P.Val < p))
  if(i == 4){
    plotMD(tfit, column=1, status=dt[, 1], main=colnames(tfit)[1], xlim=c(-8,13))
  }
  else if(i == 5){
    FirstVsSecond <- FirstVsSecond1$ENSEMBL[1:n.value]
    k <- which(v$genes$ENSEMBL %in% FirstVsSecond)
    mycol <- colorpanel(1000,"blue","white","red")
    heatmap.2(lcpm[k,], scale="row",
              labRow=v$genes$SYMBOL[k], labCol=group,
              col=mycol, trace="none", density.info="none",
              margin=c(8,8), dendrogram="column")
  }
}



ui <- fluidPage(theme = shinytheme("cerulean"),
  navbarPage("Differential Gene Expression for Placenta",



             tabPanel("Overview",
                      titlePanel(h2("Differential Gene Expression of Human Placenta Samples Using RNA Seq")),
                      br(),
                      p(
                        h4("The following analysis was done on 14 human placental samples - 8 first trimester and 6 second trimester,
                        using RNA seq data from the given dataset -", a("RNA-seq of human first and second trimester placenta.", href = "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9203")
                      )),
                      br(),
                      p(
                        h4("The allignment of samples to get counts of each gene was done using", code("STAR-seq"), "Differential gene expression
                        analysis was done using the packages -", code("edgeR"), "and", code("limma"), "from the BioConductor development software
                        for the analysis of Genomic data."
                      )),
                      br(),
                      p(h4(
                        "The webpage has different tabs for visualising", strong("Data Pre Processing, Trimester Wise Expression, DEG Analysis, Ontologies, and Pathways"),
                        "You can navigate to these tabs and select the inputs to get the relevant plots."
                      ))
                      ),


             tabPanel("Data Pre Processing",
                      sidebarPanel(
                        selectInput("type", "Type", c("MDS", "Box Plot", "Raw Data", "Filtered Data"))
                      ),
                      mainPanel(plotOutput("analysis"))
             ),

             tabPanel("Trimester Wise Expression",
                      sidebarPanel(
                        selectInput("typetrim", "Type", c("Trimester 1", "Trimester 2", "Venn Diagram")),
                        sliderInput("rangetrim", "Range", min = 0, max = 100, value = 50)
                      ),
                      mainPanel(plotOutput("plot_trim"),
                                p(h6(em(
                                  "Top expressed genes on the basis of log(CPM) values."
                                ))))
             ),

             tabPanel("DEG Analysis",
                      sidebarPanel(
                        selectInput("type1", "Type", c("Heat Map", "Diff Expressed Genes")),
                        sliderInput("p.value", "P Value", min = 0, max = 0.1, value = 0.05),
                        textOutput("slidertext")

                      ),
                      mainPanel(plotOutput("deg"))

             ),

                           tabPanel("Ontologies",
                                    sidebarPanel(
                                      selectInput("Ontology", "Ontology", c("Biological Process", "Cellular Component", "Molecular Function")),
                                      colourInput("low_color_onto", "Select low color of gradient", "black"),
                                      colourInput("high_color_onto", "Select high color of gradient", "red"),
                                      selectInput("font", "Select Font", c("Arial", "Times New Roman", "Helvetica"))
                                    ),


                                    mainPanel(plotOutput("plot_onto"),
                                              p(h6(em(
                                                "Combined score is computed by taking the log of the p-value from the Fisher exact and multiplying
                     thay by z-score of the deviation of the expected rank."
                                              ))))
                           ),



                    tabPanel("Pathways",
                     sidebarPanel(
                       selectInput("Pathway", "Pathway", c("Kegg 2021 Human", "WikiPathway 2021 Human")),
                       colourInput("low_color_path", "Select low color of gradient", "black"),
                       colourInput("high_color_path", "Select high color of gradient", "blue"),
                       selectInput("font1", "Select Font", c("Arial", "Times New Roman", "Helvetica"))
                                 ),
                                 mainPanel(plotOutput("plot_path"),
                                           p(h6(em(
                                             "Combined score is computed by taking the log of the p-value from the Fisher exact and multiplying
                                             thay by z-score of the deviation of the expected rank."
                                           )))),

     ),
     tabPanel("Genes",
              sidebarPanel(
                selectInput("gene", "select gene", "Names"),
                colourInput("first_color", "Select color of first trimester", "red"),
                colourInput("second_color", "Select color of second trimester", "blue")
              ),

              mainPanel(plotOutput("plotgenes"))),
     ))


server <- function(input, output, session) {



  output$plot_onto <- renderPlot({
    onto <- switch(input$Ontology,
                   "Biological Process" = GoPlotBio,
                   "Cellular Component" = GoPlotCell,
                   "Molecular Function" = GoPlotMol)

   plotOntology(onto, input$low_color_onto, input$high_color_onto, input$font)
    })

  output$plot_trim <- renderPlot({
     trim <- switch(input$typetrim,
                   "Trimester 1" = 1,
                   "Trimester 2" = 2,
                   "Venn Diagram" = 3)

    trimester(trim, input$rangetrim)
  })

  output$plot_path <- renderPlot({
    path <- switch(input$Pathway,
                   "Kegg 2021 Human" = kegg,
                   "WikiPathway 2021 Human" = wiki,
                   )

    plotOntology(path, input$low_color_path, input$high_color_path, input$font1)
  })


output$analysis <- renderPlot({
  n <- switch(input$type,
                 "MDS" = 1,
                 "Box Plot" = 2,
                 "Raw Data" = 3,
                 "Filtered Data" = 4
                 )
  MDSplot(n)
})

output$plotgenes <- renderPlot({
  geneplot(input$gene, input$first_color, input$second_color)
})

output$deg <- renderPlot({
  n <- switch(input$type1,
              "Diff Expressed Genes" = 4,
              "Heat Map" = 5)
  DEGgraph(n, input$p.value)
})

output$slidertext <- renderText({
  paste("The number of differentially expressed genes at p value =", input$p.value, "is", length(which(FirstVsSecond1$adj.P.Val < input$p.value)))
})

output$img <- downloadHandler(
  filename = "test.png",
  content = function(file) {
    ggsave(file, plot = DEGgraph(n, input$range), device = "png")

  })

updateSelectInput(session, inputId = "gene", choices = rownames2)

}

shinyApp(ui = ui, server = server)