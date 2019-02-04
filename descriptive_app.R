

#
uniplot <- function(x){
  if( is.factor(x) ){
    barplot(table(x)) 
  }else{
    hist(x)    
  }
}

#
biplot <- function(x, y){
  if( is.factor(x) ){
    boxplot(y~x)
  }else{
    plot(x,y)
  }
}

convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  mi
}

# ui
descriptive_app_ui <- function(dataset_names){
  
  dashboardPage(
    header = dashboardHeader(title = 'Descriptive Analysis'),
    sidebar = dashboardSidebar(
      selectInput(inputId = 'dataset_name', label = 'Dataset', choices = dataset_names),
      convertMenuItem(tabName = 'table_df', menuItem(text = 'Look@Table', tabName = 'table_df')),
      convertMenuItem(tabName = 'uniplot_df', menuItem(text = 'UniPlot', tabName = 'uniplot_df')),
      convertMenuItem(tabName = 'biplot_df', menuItem(text = 'BiPlot', tabName = 'biplot_df'))
    ),
    body = dashboardBody(
      tabItems(
        tabItem(tabName = 'table_df',
                dataTableOutput(outputId = 'head_df')),
        tabItem(tabName = 'uniplot_df',
                selectInput(inputId = 'uniplot_var', label = 'Variable', choices = ''),
                plotOutput(outputId = 'uniplot_df', width = '700px', height = '700px')),
        tabItem(tabName = 'biplot_df',
                selectInput(inputId = 'biplot_var_x', label = 'Variable X', choices = ''),
                selectInput(inputId = 'biplot_var_y', label = 'Variable Y', choices = ''),
                plotOutput(outputId = 'biplot_df', width = '700px', height = '700px'))
      )
    )
  )
  
}

# server
descriptive_app_server <- function(input, output, session){
  
  dll <- reactive({
    get(input$dataset_name)
  })
  
  output$head_df <- renderDataTable({
    df <- dll()
    df
  })
  
  get_dll_colnames <- reactive({
    colnames(dll())
  })
  
  observe({
    cn <- get_dll_colnames()
    updateSelectInput(session, inputId = 'uniplot_var', choices = cn)
    updateSelectInput(session, inputId = 'biplot_var_x', choices = cn)
  })
  
  observe({
    cn <- get_dll_colnames()
    cn <- cn[cn != input$biplot_var_x]
    updateSelectInput(session, inputId = 'biplot_var_y', choices = cn)
  })
  
  output$uniplot_df <- renderPlot({
    uniplot(dll()[,input$uniplot_var])
  })
  
  output$biplot_df <- renderPlot({
    biplot(x = dll()[,input$biplot_var_x], y = dll()[,input$biplot_var_y])
  })
  
}

# launch ShinyApp
descriptive_app <- function(){
  
  # get dataframe candidates
  workspace_elem_list <- ls(envir = .GlobalEnv)
  elem_is_dataframe <- sapply( workspace_elem_list, function(e){
    is.data.frame(get(e))
  })
  cand <- workspace_elem_list[elem_is_dataframe]
  
  # launch app
  shinyApp(ui = descriptive_app_ui(dataset_names = cand), server = descriptive_app_server)
}
descriptive_app()

