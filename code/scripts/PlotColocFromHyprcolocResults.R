#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : ColocPlots
# @created     : Thursday May 26, 2022 10:49:00 CDT
#
# @description : Script to make coloc plots from my custom hyprcoloc tsv
# output. Run in code dir. some input files are hardcoded
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    # args <- scan(text=
    #              "scratch/TestHyprcolocPlots.tsv scratch/ColocPlotTest_", what='character')
    args <- scan(text="scratch/test_hyprcoloc.txt.gz scratch/ColocExamples/Group_MetabolicAfterClumpingSnps.Plot. png", what='character')
    # args <- scan(text="scratch/ForPoster/PotentialColocExamples.tsv scratch/ForPoster/Coloc_ png", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

hyprcoloc_output_f_in <- args[1]
figure_out_prefix <- args[2]
figure_out_type <- args[3] #either pdf or png


library(tidyverse)
library(data.table)
library(GGally)
library(ggrepel)
library(RColorBrewer)

color_pal_fun <- colorRampPalette(brewer.pal(7,"Dark2"))

lower_point <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(..., size=0.5, alpha=0.2, color="gray") + 
    geom_point(data = (
                       data %>% filter(!is.na(TopCandidateSNP))),
                       aes(fill=Color),
                       shape = 21, color="black", alpha=1, size=1.5
                       ) +
    scale_fill_identity() +
    # scale_fill_brewer(palette="Dark2", type="qualitative", na.value="gray") +
    theme_classic()
}
my_fn <- function(data, mapping, method="p", use="pairwise", ...){
              # grab data
              x <- eval_data_col(data, mapping$x)
              y <- eval_data_col(data, mapping$y)
              # calculate correlation
              corr <- cor(x, y, method=method, use=use)
              # calculate colour based on correlation value
              # Here I have set a correlation of minus one to blue, 
              # zero to white, and one to red 
              # Change this to suit: possibly extend to add as an argument of `my_fn`
              colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
              fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
              ggally_cor(data = data, mapping = mapping, ...) + 
                theme_void() +
                theme(panel.background = element_rect(fill=fill))
            }

ExpressedGenes <- fread("QTLs/QTLTools/Expression.Splicing/PermutationPassForColoc.txt.gz", select=1:4) %>%
    mutate(phe_id = str_replace(phe_id, "(.+?):.+?$", "\\1"))

# dat %>% mutate(Color = recode(TopCandidateSNP, !!!(unique(dat$TopCandidateSNP) %>% setNames( color_pal_fun(length(.)), . )), .default="gray"))

AllPhenotypesTrueLocation <- Sys.glob("QTLs/QTLTools/*/OnlyFirstReps.sorted.qqnorm.bed.gz") %>%
    setNames(str_replace(., "QTLs/QTLTools/(.+?)/OnlyFirstReps.sorted.qqnorm.bed.gz", "\\1")) %>%
    lapply(fread, select=1:6) %>%
    bind_rows(.id = "phe_class") %>%
    unite(Trait, phe_class, pid, sep=";") %>%
    mutate(x = case_when(
                         strand == "-" ~ end,
                         TRUE ~ start
                         )) %>%
    mutate(xend = case_when(
                            strand == "-" ~ start,
                            TRUE ~ end
                            ))

dat.all <- fread(hyprcoloc_output_f_in) %>%
    inner_join(ExpressedGenes, by=c("GeneLocus"="phe_id")) %>%
    separate(Trait, into = c("phe_class","phe_id"), sep=";") %>%
    base::split(.$GeneLocus)


for (dat in dat.all){
    ColorsMap <- dat %>% pull(TopCandidateSNP) %>% unique() %>%
        na.omit() %>% sort() %>%
        setNames(color_pal_fun(length(.)), .)
    summarystats <- dat %>%
        distinct(GeneLocus, phe_class,  .keep_all=T) %>%
        mutate(cmd = str_glue("tabix -h QTLs/QTLTools/{phe_class}/NominalPassForColoc.txt.tabix.gz {phe_chr}:{phe_from}-{phe_to}")) %>%
        select(phe_class, cmd) %>%
        deframe() %>%
        lapply(fread, sep='\t', select=c(1, 8, 10, 12), colClasses = c("#phe_id"="character", "var_id"="character", "var_from"="numeric", "nom_pval"="numeric")) %>%
        bind_rows(.id = "phe_class") %>%
        mutate(phe_id_GeneLocus = str_replace(`#phe_id`, "(.+):(.+?)$", "\\1;\\2")) %>%
        separate(phe_id_GeneLocus, into = c("phe_id", "GeneLocus"), sep=";") %>%
        select(GeneLocus, var_from, nom_pval, phe_id, phe_class, var_id) %>%
        inner_join(dat, by=c("GeneLocus", "phe_id", "phe_class"))
    ToPlot <- summarystats %>%
        unite(Trait, phe_class, phe_id, sep=";") %>%
        # rename traits here
        mutate(Sig = -log10(nom_pval)) %>%
        select(Trait, Sig, var_id, var_from, TopCandidateSNP) %>%
        arrange(TopCandidateSNP) %>%
        # distinct(.keep_all=T) %>%
        pivot_wider(names_from = Trait, values_from = Sig, id_cols = c("var_id", "var_from")) %>%
        drop_na() %>%
        mutate(TopCandidateSNP = case_when(
                                           var_id %in% dat$TopCandidateSNP ~ var_id,
                                           TRUE ~ NA_character_
                                             )) %>%
        mutate(Color = recode(TopCandidateSNP, !!!ColorsMap, .default="gray")) %>%
        select(var_id, var_from, TopCandidateSNP, Color, everything()) %>%
        arrange(desc(TopCandidateSNP))

    P <- ggpairs(ToPlot,
            columnLabels = gsub('[;.]', ' ', colnames(ToPlot)[-c(1:4)], perl=T), 
            labeller = label_wrap_gen(10),
            columns = 5:ncol(ToPlot),
            diag=list(continuous = "blankDiag"),
            lower=list(continuous = lower_point),
            upper=list(continuous = my_fn)
            )
    SymetricMatrixLabels <- colnames(ToPlot)[-c(1:4)]
    for (i in seq_along(SymetricMatrixLabels)){
        RectangleDat <- dat %>%
            unite(Trait, phe_class, phe_id, sep=";") %>%
            filter(Trait == SymetricMatrixLabels[i]) %>%
            mutate(Color = recode(TopCandidateSNP, !!!ColorsMap, .default="gray"))
        SegmentDat_unstranded <- AllPhenotypesTrueLocation %>%
            filter(Trait == SymetricMatrixLabels[i]) %>%
            filter(!str_detect(strand, "[+-]"))
        SegmentDat_stranded <- AllPhenotypesTrueLocation %>%
            filter(Trait == SymetricMatrixLabels[i]) %>%
            filter(str_detect(strand, "[+-]"))
        # To continuehere 
        ToPlotForManhattan <- ToPlot %>%
            select(var_from, SymetricMatrixLabels[i], TopCandidateSNP, Color) %>%
            rename("Sig"=2)
        P[i,i] <- ggplot(ToPlotForManhattan) +
            geom_rect(data = RectangleDat, aes(color=Color), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0, size=2) +
            geom_point(aes(x=var_from/1E6, y=Sig), size=0.5, alpha=0.2, color="gray") +
            geom_segment(data = SegmentDat_stranded, aes(x=x/1E6, xend=xend/1E6), y=-0.1, yend=-0.1, arrow = arrow(length = unit(0.1,"cm"))) +
            geom_segment(data = SegmentDat_unstranded, aes(x=x/1E6, xend=xend/1E6), y=-0.1, yend=-0.1) +
            geom_point(data = (
                       ToPlotForManhattan %>% filter(!is.na(TopCandidateSNP))),
                       aes(x=var_from/1E6, y=Sig, fill=Color),
                       shape = 21, color="black", alpha=1, size=1.5
                       ) +
            scale_fill_identity() +
            scale_color_identity() +
            # scale_fill_brewer(palette="Dark2", type="qualitative", na.value="gray") +
            theme_bw() +
            theme(
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            )
            if (i==length(SymetricMatrixLabels)){
                P[i,i] <- P[i,i] +
                      coord_cartesian(clip="off") +
                      geom_text(label = paste("Position (Mb)", RectangleDat$phe_chr), hjust=0, vjust=4, y=-Inf, x=RectangleDat$phe_from/1E6, size=2) +
                      theme(plot.margin = margin(b=25))
            }
    }
    P_fn <- paste0(figure_out_prefix, dat$GeneLocus[1], ".", figure_out_type)
    print(str_glue("saving plot file: {P_fn}"))
    ggsave(P_fn, P, height=length(SymetricMatrixLabels)*1.2, width=length(SymetricMatrixLabels)*1.2, limitsize=F)
}

