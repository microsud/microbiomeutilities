
plot_spaghetti_taxa <- function(x,
                                select.taxa=NULL,
                                xvar=NULL,
                                group,
                                line.bg.color="#8d99ae",
                                bg.opacity=0.5,
                                focus.color="brown3",
                                ncol=NULL,
                                nrow=NULL,
                                focus.line.size=0.5,
                                line.size=1) {
  
  if (is.null(xvar)){
    stop("xvar cannot be empty")
  }
  OTUID <- name2 <- NULL
  phdf <- phy_to_ldf(x, transform.counts = NULL)
  
  # https://www.data-to-viz.com/caveat/spaghetti.html
  tmp <- phdf %>%
    mutate(name2=OTUID)
  
  tmp %>%
    ggplot(aes_string(x=xvar, y="Abundance")) +
    geom_line(data=tmp %>% dplyr::select(-OTUID), 
              aes(group=name2), 
              color=line.bg.color, size=line.size, alpha=bg.opacity) +
    geom_line(aes(color=OTUID), color=focus.color, size=focus.line.size)+
    theme_biome_utils() +
    theme(
      legend.position="none"
    ) +
    facet_wrap(~OTUID, 
               ncol = ncol,
               nrow=nrow)
}


################################################################################


plot_spaghetti_sample <- function(x,
                                  select.taxa=NULL,
                                  xvar=NULL,
                                  group,
                                  line.bg.color="#8d99ae",
                                  bg.opacity=0.5,
                                  focus.color="brown3",
                                  ncol=NULL,
                                  nrow=NULL,
                                  focus.line.size=0.5,
                                  line.size=1){
  
  if (is.null(xvar)){
    stop("xvar cannot be empty")
  }
  
  OTUID <- group2 <- name2 <- select.tax <- NULL
  
  phdf <- phy_to_ldf(x, transform.counts =NULL) %>% 
    filter(OTUID==select.taxa)
  
  # https://www.data-to-viz.com/caveat/spaghetti.html
  sub.var <- sym(group)
  
  tmp <- phdf %>%
    mutate(name2=OTUID,
           group2=!!sub.var)
  
  tmp %>%
    ggplot(aes_string(x=xvar, y="Abundance")) +
    geom_line(data=tmp %>% dplyr::select(-group2), 
              aes(group=!!sub.var), 
              color=line.bg.color, size=line.size, alpha=bg.opacity) +
    geom_line(aes(color=OTUID), color=focus.color, size=focus.line.size)+
    theme_biome_utils() +
    theme(
      legend.position="none"
    ) +
    facet_wrap(~group2, 
               ncol = ncol,
               nrow=nrow)
  
}

