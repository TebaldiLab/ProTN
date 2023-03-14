# ----------------- MAIN REACTIVE FUNCTION ---------------- #

newdf <- function(comp, kinase_df) {
  
  # get current values
  tempdf = svginfo$dataframe
  
  # set font family
  tempdf$text.font = "'Helvetica'"
  
  # establish legend
  legend = c()
  # Set initial yoffset
  yoffset = 79.125
  
  # get current values
  tempdf$text.size = 4
  
  {
    # read in text area input
    recolordf = data.frame(
      "kinase" = names(kinase_df),
      "userinfo" = kinase_df,
      stringsAsFactors = F
    )
    
    # convert to coral id
    recolordf = convertID (tempdf, recolordf, inputtype = "HGNC")
    
    if (nrow(recolordf) > 0)
    {
      {
        branchcolorpalette_div <- "Blue_Grey_Coral"
        palettecolors = unlist(divpalettes[branchcolorpalette_div])
        branchcolpalette = colorRampPalette(palettecolors)(11)
        
        bg.col = unlist(divpalettes[branchcolorpalette_div])[2]
      }
      
      # set colors based on group
      newcolors_and_colormapping = color.by.value(
        df = tempdf,
        recolordf = recolordf,
        colors  = branchcolpalette,
        heatrange = c(as.integer(min(
          recolordf$values
        )) - 1, as.integer(max(
          recolordf$values
        )) + 1),
        bg.col = bg.col
      )
      tempdf$branch.col = newcolors_and_colormapping[[1]]
      tempdf$branch.val = newcolors_and_colormapping[[2]]
      
      # reorder based on branch value
      tempdf$branchorder = order(abs(tempdf$branch.val),
                                 decreasing = FALSE,
                                 na.last = FALSE)
      
      # add legend info
      lines_and_offset = build.value.legend(
        yoffset = yoffset,
        minval = as.integer(min(recolordf$values)) - 1,
        maxval = as.integer(max(recolordf$values)) + 1,
        palette = branchcolpalette,
        elementtype = "Branch",
        fontfamily = tempdf$text.font,
        subtitle = "Kinase Activity Score"
      )
      lines = lines_and_offset[[1]]
      yoffset = lines_and_offset[[2]] + 14
      legend = c(legend, lines)
    }
  }
  
  # ------------------ NODE COLOR ------------------ #
  {
    # read in text area input
    recolordf = data.frame(
      "kinase" = names(kinase_df),
      "userinfo" = kinase_df,
      stringsAsFactors = F
    )
    
    # convert to coral id
    recolordf = convertID (tempdf, recolordf, inputtype = "HGNC")
    
    if (nrow(recolordf) > 0)
    {
      {
        nodecolorpalette_div <- "Blue_Grey_Coral"
        palettecolors = unlist(divpalettes[nodecolorpalette_div])
        nodecolpalette = colorRampPalette(palettecolors)(11)
        bg.col = unlist(divpalettes[nodecolorpalette_div])[2]
      }
      {
        bg.col = "#FFFFFF"
      }
      
      # set colors based on value
      newcolors_and_colormapping = color.by.value(
        df = tempdf,
        recolordf = recolordf,
        colors  = nodecolpalette,
        heatrange = c(as.integer(min(
          recolordf$values
        )) - 1, as.integer(max(
          recolordf$values
        )) + 1),
        bg.col = bg.col
      )
      tempdf$node.col = newcolors_and_colormapping[[1]]
      tempdf$node.val = newcolors_and_colormapping[[2]]
      
      # reorder based on branch color
      tempdf$nodeorder = order(abs(tempdf$node.val),
                               decreasing = FALSE,
                               na.last = FALSE)
    }
  }
  
  # ------------------ NODE SIZE ------------------ #
  
  # color nodes by single color
  # if (input$nodesizetype == "One Size")
  {
    tempdf$node.radius = 5
  }
  # ------------------ ADVANCED OPTIONS ------------------ #
  
  tempdf$node.opacity = 1
  
  # Change Color and Size of Font for Selected kinases
  # if (input$fontcolorselect == "Manual")
  {
    selkinases = names(kinase_df)
    
    selkinasescoral = ""
    if (length(selkinases) > 0)
    {
      # convert selected to coral ids
      kinasestoconvert = data.frame(kin1 = selkinases, kin2 = selkinases)
      selkinasesconverted = convertID (tempdf, kinasestoconvert, inputtype =
                                         "HGNC")
      if (nrow(selkinasesconverted) > 0)
      {
        selkinasescoral = selkinasesconverted[, 1]
      }
    }
    
    # set background color and font size
    tempdf$text.col = "#FFFFFF"
    tempdf$text.size = 0.15
    
    # set selected font color and size
    tempdf$text.col[which(tempdf$id.coral %in% selkinasescoral)]  = "#000000"
    tempdf$text.size[which(tempdf$id.coral %in% selkinasescoral)] = 6
  }
  
  # if (input$nodestrokecolselect == "Same as Node")
  {
    tempdf$node.strokecol = tempdf$node.col
  }
  
  return(list(tempdf, legend))
} # end reactive

# ----------------- PLOTS ---------------- #

# build the manning tree
renderSvgPanZoom <- function(comp, kinase_Act, dirOutput_kinase) {
  # recolor the official matrix
  dfandlegend = newdf(comp, kinase_Act)
  svginfo$dataframe = dfandlegend[[1]]
  svginfo$legend = dfandlegend[[2]]
  
  # set title
  svginfo$title = paste0("Kinase Tree - Comparison: ", comp)
  
  # Write SVG file
  svgoutfile <- writekinasetree(
    svginfo,
    destination = paste0(dirOutput_kinase, comp, "_kinase_Tree_CORAL.svg"),
    font = svginfo$dataframe$text.font[1],
    labelselect = "HGNC",
    groupcolor = "#000000"
  )
  # Render SVG
  return(paste0(dirOutput_kinase, comp, "_kinase_Tree_CORAL.svg"))
  # svgPanZoom(svgoutfile,viewBox = T,controlIconsEnabled=F)
}


# ----------------- DATA TABLE ---------------- #

# build the table
# output$KinaseTable <- DT::renderDataTable({
#
#   dfandlegend = newdf()
#   simpldf = dfandlegend[[1]][,c("id.coral","id.longname","kinase.family","kinase.group","branch.val","branch.col","node.col","node.radius")]
#
#   # reverse the data frame so that colored kinases appear first
#   simpldf<-simpldf[dim(simpldf)[1]:1,]
#
#   # convert branch colors to rgb
#   mycolors <- simpldf$branch.col
#   rgbcolors <- apply(grDevices::col2rgb(mycolors), 2,
#                      function(rgb) sprintf("rgb(%s)", paste(rgb, collapse=",")))
#   tgtbranch <- sprintf('<span style="color:%s">&#9608;</span>', rgbcolors)
#
#   newdf = data.frame(kinase=simpldf$id.coral,
#                      name=simpldf$id.longname,
#                      family=simpldf$kinase.family,
#                      group = simpldf$kinase.group
#   )
#
#   # add node info if present
#   if ("none" %in% simpldf$node.col == F)
#   {
#     # convert node colors to rgb
#     mycolors <- simpldf$node.col
#
#     rgbcolors <- apply(grDevices::col2rgb(mycolors), 2,
#                        function(rgb) sprintf("rgb(%s)", paste(rgb, collapse=",")))
#     tgtnode <- sprintf('<span style="color:%s">&#11044;</span>', rgbcolors)
#
#     newdf$node_radius = simpldf$node.radius
#     newdf$node_color = tgtnode
#   }
#
#   # add branch color
#   newdf$branch_color=tgtbranch
#
#   datatable(newdf, escape=FALSE)
# })
#
