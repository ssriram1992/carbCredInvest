#Leer las librerias requeridas
library(MASS)
library(ISLR)
library(ggplot2)
library(caret)
library(tidyr)
library(psych)
library(corrplot)
library(plotly)
###############################
library(readxl)
library(jtools)
library(sjPlot)
library(sjmisc)
library(effects)
library(gridExtra)
library(grid)
library(ggiraphExtra)
library(ggeffects)
library(ggthemes)
library(ggpubr)

##########################################################################################
data <- read_excel("MPECwithDecrCarbCred_data.xlsx")
data$Trade<-"Trade"
dataNt <- read_excel("MPECwithDecrCarbCredNoTrade_data.xlsx")
dataNt$Trade<-"NoTrade"



dataEPEC <- read_excel("EPEC_analysis_Data2.xlsx")
dataEPEC$Trade<-"EPEC"

plot_g<-function(y,x,lab_y,lab_x,y_low,y_high){
  ggplot(data, aes(y=y,x=x, fill=Trade))+
    geom_boxplot(outlier.shape = NA)+
    #theme_set(theme_sjplot())+ 
    theme_set(theme_bw())+
    theme(plot.title=element_blank(),
          axis.text=element_text(size=12),
          legend.position="bottom")+
    xlab(lab_x)+
    ylab(lab_y)+
    ylim(y_low,y_high)+ 
    scale_fill_brewer(palette="Pastel2")
}

plot_facet<-function(y,x,lab_y,lab_x,y_low,y_high){
  ggplot(data, aes(y=y,x=x, fill=Trade))+
    geom_boxplot(outlier.shape = NA)+
    #theme_set(theme_sjplot())+
    theme_set(theme_bw())+
    theme(plot.title=element_blank(),
          axis.text=element_text(size=12),
          #legend.position="bottom",
          legend.position = "none")+
    xlab(lab_x)+
    ylab(lab_y)+
    ylim(y_low,y_high)+ 
    scale_fill_brewer(palette="Pastel2")+
    facet_wrap(.~Trade)
}

plot_facet_e<-function(y,x,lab_y,lab_x,y_low,y_high){
  ggplot(dataEPEC, aes(y=y,x=x, fill=Trade))+
    geom_boxplot(outlier.shape = NA)+
    #theme_set(theme_sjplot())+
    theme_set(theme_bw())+
    theme(plot.title=element_blank(),
          axis.text=element_text(size=12),
          #legend.position="bottom",
          legend.position = "none")+
    xlab(lab_x)+
    ylab(lab_y)+
    ylim(y_low,y_high)+ 
    scale_fill_brewer(palette="Pastel2")+
    facet_wrap(.~Trade)
}

data<-rbind(data,dataNt)
names(data)
data$Trade<-as.factor(data$Trade)
data$inVS<-as.factor(data$inVS)
data$demS<-as.factor(data$demS)
data$sdCF<-as.factor(data$sdCF)
data$eI_c1<-as.factor(data$eI_c1)
data$eI_c2<-as.factor(data$eI_c2)
data$crc1<-as.factor(data$crc1)
data$crc2<-as.factor(data$crc2)
data$init_A<-as.factor(data$init_A)
names(data)
#####################################################################################################

a<-plot_g(data$CPriceC1,data$inVS," (MWh)", "inVS",0,500)
b<-plot_g(data$CPriceC1,data$demS," (MWh)", "Demand price sensitivity",0,500)
c<-plot_g(data$CPriceC1,data$sdCF," (MWh)", "St. Dev. Capacity Factor",0,500)
d<-plot_g(data$CPriceC1,data$eI_c1," (MWh)", "eI_c1",0,500)
e<-plot_g(data$CPriceC1,data$eI_c2," (MWh)", "eI_c2",0,500)
f<-plot_g(data$CPriceC1,data$crc1," (MWh)", "crc1 (same as crc2)",0,500)
g<-plot_g(data$CPriceC1,data$crc2," (MWh)", "crc1",0,500)


a<-plot_facet(data$TotalEnergy,data$inVS,"Total Energy (MWh)", "inVS",0,500)
b<-plot_facet(data$TotalEnergy,data$demS,"Total Energy (MWh)", "Demand price sensitivity",0,500)
c<-plot_facet(data$TotalEnergy,data$sdCF,"Total Energy (MWh)", "St. Dev. Capacity Factor",0,500)
d<-plot_facet(data$TotalEnergy,data$eI_c1,"Total Energy (MWh)", "eI_c1",0,500)
e<-plot_facet(data$TotalEnergy,data$eI_c2,"Total Energy (MWh)", "eI_c2",0,500)
f<-plot_facet(data$TotalEnergy,data$crc1,"Total Energy (MWh)", "crc1 (same as crc2)",0,500)
g<-plot_facet(data$TotalEnergy,data$crc2,"Total Energy (MWh)", "crc1",0,500)



#1200x740
grid.arrange(a, b, c, d, e, f,ncol=3+ theme(legend.position="none"))
grid.arrange(a, b, c, d, e, f,ncol=3+ theme(legend.position="none"))
ggarrange(a, b, c, d, e, f,ncol=3, nrow=2,common.legend = TRUE,legend="bottom")

annotate_figure(ggarrange(a, b, c, d, e, f,ncol=3, nrow=2,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total Energy Generation", size = 16))
######################################################################################################


a<-plot_g(data$TotalRenew,data$inVS,"(MWh)", "inVS",100,250)
b<-plot_g(data$TotalRenew,data$demS,"(MWh)", "Demand price sensitivity",100,250)
c<-plot_g(data$TotalRenew,data$sdCF,"(MWh)", "St. Dev. Capacity Factor",100,250)
d<-plot_g(data$TotalRenew,data$eI_c1,"(MWh)", "eI_c1",100,250)
e<-plot_g(data$TotalRenew,data$eI_c2,"(MWh)", "eI_c2",100,250)
f<-plot_g(data$TotalRenew,data$crc1,"(MWh)", "crc1 (same as crc2)",100,250)
g<-plot_g(data$TotalRenew,data$crc2,"(MWh)", "crc1",100,250)

a<-plot_facet(data$TotalRenew,data$inVS,"Total green Energy (MWh)", "inVS",100,250)
b<-plot_facet(data$TotalRenew,data$demS,"Total green Energy (MWh)", "Demand price sensitivity",100,250)
c<-plot_facet(data$TotalRenew,data$sdCF,"Total green Energy (MWh)", "St. Dev. Capacity Factor",100,250)
d<-plot_facet(data$TotalRenew,data$eI_c1,"Total green Energy (MWh)", "eI_c1",100,250)
e<-plot_facet(data$TotalRenew,data$eI_c2,"Total green Energy (MWh)", "eI_c2",100,250)
f<-plot_facet(data$TotalRenew,data$crc1,"Total green Energy (MWh)", "crc1 (same as crc2)",100,250)
g<-plot_facet(data$TotalRenew,data$crc2,"Total green Energy (MWh)", "crc1",100,250)

grid.arrange(a, b, c, d, e, f,ncol=3)
#1200x735
annotate_figure(ggarrange(a, b, c, d, e, f,ncol=3, nrow=2,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total Green Generation", size = 16))

#######################################################################################################
a<-plot_g(data$TotalDirty,data$inVS,"(MWh)", "inVS",0,300)
b<-plot_g(data$TotalDirty,data$demS,"(MWh)", "Demand price sensitivity",0,300)
c<-plot_g(data$TotalDirty,data$sdCF,"(MWh)", "St. Dev. Capacity Factor",0,300)
d<-plot_g(data$TotalDirty,data$eI_c1,"(MWh)", "eI_c1",0,300)
e<-plot_g(data$TotalDirty,data$eI_c2,"(MWh)", "eI_c2",0,300)
f<-plot_g(data$TotalDirty,data$crc1,"(MWh)", "crc1 (same as crc2)",0,300)
g<-plot_g(data$TotalDirty,data$crc2,"(MWh)", "crc1",0,300)


grid.arrange(a, b, c, d, e, f,ncol=3)
#1200x735

a<-plot_facet(data$TotalDirty,data$inVS,"(MWh)", "inVS",0,300)
b<-plot_facet(data$TotalDirty,data$demS,"(MWh)", "Demand price sensitivity",0,300)
c<-plot_facet(data$TotalDirty,data$sdCF,"(MWh)", "St. Dev. Capacity Factor",0,300)
d<-plot_facet(data$TotalDirty,data$eI_c1,"(MWh)", "eI_c1",0,300)
e<-plot_facet(data$TotalDirty,data$eI_c2,"(MWh)", "eI_c2",0,300)
f<-plot_facet(data$TotalDirty,data$crc1,"(MWh)", "crc1 (same as crc2)",0,300)
g<-plot_facet(data$TotalDirty,data$crc2,"(MWh)", "crc1",0,300)


grid.arrange(a, b, c, d, e, f,
             ncol=3,
             top = textGrob("Total Fossil Energy",
                            gp=gpar(fontsize=16)))
#1200x740
annotate_figure(ggarrange(a, b, c, d, e, f,ncol=3, nrow=2,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total Fossil Generation", size = 16))
####################################################
##########EPEC

str(dataEPEC)
names(dataEPEC)

dataEPEC$suppInt<-as.factor(dataEPEC$suppInt)
dataEPEC$suppSlope<-as.factor(dataEPEC$suppSlope)
dataEPEC$c1_prodnVal<-as.factor(dataEPEC$c1_prodnVal)
dataEPEC$c1_Tax<-as.factor(dataEPEC$c1_Tax)
dataEPEC$c2_Tax<-as.factor(dataEPEC$c2_Tax)
dataEPEC$c2_prodnVal<-as.factor(dataEPEC$c2_prodnVal)

#dataEPEC<-dataEPEC[dataEPEC$TotalProd<=1000,]
#dataEPEC<-dataEPEC[dataEPEC$TotalProd>=500,]
plot(dataEPEC$TotalProd)
plot(dataEPEC$Supply)
plot(dataEPEC$Price)
plot(dataEPEC$Emission)
plot(dataEPEC$SolInv)
plot(dataEPEC$WindInv)
plot(dataEPEC$dirtyProd)
plot(dataEPEC$CleanProd)

names(dataEPEC)
a<-plot_facet_e(dataEPEC$TotalProd,dataEPEC$suppInt,"(MWh)", "suppInt",0,1000)
b<-plot_facet_e(dataEPEC$TotalProd,dataEPEC$suppSlope,"(MWh)", "suppSlope",0,1000)
c<-plot_facet_e(dataEPEC$TotalProd,dataEPEC$c1_Tax,"(MWh)", "c1_Tax (same for c2_Tax)",0,1000)
d<-plot_facet_e(dataEPEC$TotalProd,dataEPEC$c2_Tax,"(MWh)", "c2_Tax",0,1000)
e<-plot_facet_e(dataEPEC$TotalProd,dataEPEC$c1_prodnVal,"(MWh)", "c1_prodnVal (same for c2_prodnVal)",0,1000)
f<-plot_facet_e(dataEPEC$TotalProd,dataEPEC$c2_prodnVal,"(MWh)", "c2_prodnVal",0,1000)

annotate_figure(ggarrange(a, b, c, e,ncol=4, nrow=1,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total Generation", size = 16))

#1400x400
######################################################################
a<-plot_facet_e(dataEPEC$CleanProd,dataEPEC$suppInt,"(MWh)", "suppInt",100,700)
b<-plot_facet_e(dataEPEC$CleanProd,dataEPEC$suppSlope,"(MWh)", "suppSlope",100,700)
c<-plot_facet_e(dataEPEC$CleanProd,dataEPEC$c1_Tax,"(MWh)", "c1_Tax (same for c2_Tax)",100,700)
d<-plot_facet_e(dataEPEC$CleanProd,dataEPEC$c2_Tax,"(MWh)", "c2_Tax",100,700)
e<-plot_facet_e(dataEPEC$CleanProd,dataEPEC$c1_prodnVal,"(MWh)", "c1_prodnVal (same for c2_prodnVal)",100,700)
f<-plot_facet_e(dataEPEC$CleanProd,dataEPEC$c2_prodnVal,"(MWh)", "c2_prodnVal",100,700)

annotate_figure(ggarrange(a, b, c, e,ncol=4, nrow=1,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total Green Generation", size = 16))

##1400x400

######################################################################
a<-plot_facet_e(dataEPEC$Supply,dataEPEC$suppInt,"(TCO2)", "suppInt",0,1000)
b<-plot_facet_e(dataEPEC$Supply,dataEPEC$suppSlope,"(TCO2)", "suppSlope",0,1000)
c<-plot_facet_e(dataEPEC$Supply,dataEPEC$c1_Tax,"(TCO2)", "c1_Tax (same for c2_Tax)",0,1000)
d<-plot_facet_e(dataEPEC$Supply,dataEPEC$c2_Tax,"(TCO2)",  "c2_Tax",0,1000)
e<-plot_facet_e(dataEPEC$Supply,dataEPEC$c1_prodnVal,"(TCO2)", "c1_prodnVal (same for c2_prodnVal)",0,1000)
f<-plot_facet_e(dataEPEC$Supply,dataEPEC$c2_prodnVal,"(TCO2)", "c2_prodnVal",0,1000)

annotate_figure(ggarrange(a, b, c, e,ncol=4, nrow=1,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total allowance supply", size = 16))

#plot(dataEPEC$c1_Tax,dataEPEC$TotalRenew)
######################################################################
a<-plot_facet_e(dataEPEC$Price,dataEPEC$suppInt,"($/TCO2)", "suppInt",0,400)
b<-plot_facet_e(dataEPEC$Price,dataEPEC$suppSlope,"($/TCO2)", "suppSlope",0,400)
c<-plot_facet_e(dataEPEC$Price,dataEPEC$c1_Tax,"($/TCO2)", "c1_Tax (same for c2_Tax)",0,400)
d<-plot_facet_e(dataEPEC$Price,dataEPEC$c2_Tax,"($/TCO2)", "c2_Tax",0,400)
e<-plot_facet_e(dataEPEC$Price,dataEPEC$c1_prodnVal,"($/TCO2)", "c1_prodnVal (same for c2_prodnVal)",0,400)
f<-plot_facet_e(dataEPEC$Price,dataEPEC$c2_prodnVal,"($/TCO2)", "c2_prodnVal",0,400)

annotate_figure(ggarrange(a, b, c, e,ncol=4, nrow=1,common.legend = TRUE,legend="bottom"),
                top = text_grob("Allowance price", size = 16))

#plot(dataEPEC$c1_Tax,dataEPEC$TotalRenew)
######################################################################

a<-plot_facet_e(dataEPEC$Emission,dataEPEC$suppInt,"(TCO2)","suppInt",100,700)
b<-plot_facet_e(dataEPEC$Emission,dataEPEC$suppSlope,"(TCO2)", "suppSlope",100,700)
c<-plot_facet_e(dataEPEC$Emission,dataEPEC$c1_Tax,"(TCO2)", "c1_Tax (same for c2_Tax)",100,700)
d<-plot_facet_e(dataEPEC$Emission,dataEPEC$c2_Tax,"(TCO2)", "c2_Tax",100,700)
e<-plot_facet_e(dataEPEC$Emission,dataEPEC$c1_prodnVal,"(TCO2)", "c1_prodnVal (same for c2_prodnVal)",100,700)
f<-plot_facet_e(dataEPEC$Emission,dataEPEC$c2_prodnVal,"(TCO2)", "c2_prodnVal",100,700)


annotate_figure(ggarrange(a, b, c, e,ncol=4, nrow=1,common.legend = TRUE,legend="bottom"),
                top = text_grob("Total Emissions", size = 16))
