#data preprocessing from lab absorbance reader
library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)

install.packages("devtools")
yes
devtools::install_github("JIC-CSB/babar")
library(babar)
#file necessary.
#metadata for each well spearated by columns named []_metadata.csv
#raw data files unprocessed



#Save constants of molar extinction coefficients of PrestoBlue reagent 
#these values are from manufacturer you will multiply all wells by this number

Reduction_570nm<-155677
Oxidized_570nm<-80586

Reduction_600nm<-14652
Oxidized_600nm<-117216
#import data
Meta <- read.csv("RawData/1201MM_metadata.csv")

RawData570 <-
  read.csv("RawData/1201MM_12082022113915.csv",
           skip = 43,
           nrows = 96)

RawData600<-  read.csv("RawData/1201MM_12082022113956.csv",
                       skip = 43,
                       nrows = 96)


#change time to  a time recognized character string




#change readings to time
RawData600[, 2] <- hms(RawData600[, 2])
#check one well to make sure you get a graph
plot(RawData600$Time, RawData600$A1)

#change readings to time
RawData570[, 2] <- hms(RawData570[, 2])
#check one well to make sure you get a graph
plot(RawData570$Time, RawData570$A1)

#make new data set 


#Perwell equation  = ((Oxidized_600nm*Test_Absorbance570)-(Oxidized_570nm*Test_Absorbance600))/((Reduction_570nm*Control_Absorbance600)-(Reduction_600nm*Control_Absorbance570))


i<-1

#get Control wells(Media + reagent)

ControlWells<-c("C1",'C2','C3','C4','C5','C6')

for (i in 1:nrow(RawData570)) {
  
#calculate the average for that point in time of control wells
  #we do per time point incase media causes reduction
RowStandard<-mean(unlist((RawData600[i,ControlWells]*Reduction_570nm)-(RawData570[i,ControlWells]*Reduction_600nm)))

#per Well j

for (j in 3:ncol(RawData570)){

RawData570[i,j]<-((Oxidized_600nm*RawData570[i,j])-(Oxidized_570nm*RawData600[i,j]))/RowStandard
}
}

#rename Rawdata570

Rawdata_PercentReduction<-RawData570
#check one well to make sure you get a graph
plot(Rawdata_PercentReduction$Time, Rawdata_PercentReduction$A2)




View(Rawdata_PercentReduction[1:5,1:5])



#make a long version of the dataframe
Data <- Rawdata_PercentReduction %>% pivot_longer(cols = c(Meta$Well),
                                 names_to = 'Well', values_to = 'Abs')

#merge by well ID
Data <- merge(Data, Meta, by = "Well")
Data$Well <- as.factor(Data$Well)
Data$Isolate <- as.factor(Data$Isolate)
Data$Dosage_uM<-as.factor(Data$Dosage_uM)
Data$Spores <- as.factor(Data$Spores)


#remove outiler wells
#Data<-Data[Data$Well!="E4"]
#Data <- subset(Data,!(Well %in% c("C1", "B7", "C5","H9")))

#remove Contaminated wells & blanks
#Data <- subset(Data,!(Isolate %in% c("Contam","blank")))

#Data <- subset(Data,(SporeNumber %in% c("300 spores")))
#Data <- subset(Data,(Well %in% c("D7","D8","D9")))

#remove C & F rows for plotting

CF_remove<-c(paste0("C",seq(from=1, to=12, by=1)), paste0("F",seq(from=1, to=12, by=1)))

Data <- subset(Data,!(Well %in% CF_remove))

#remove Contols
Data[1:5,1:5]

View(Data[1:5,1:5])

# multiple line plot
#convert Data to oxygen used per well 
#wave length

#we meed to find the point of absorbance = .2 
find_intersect <- function(x,y, target=0.5) {
  optimize(function(z) (approxfun(x,y)(z)-target)^2, x)$minimum
}
#segment data based on categories to be graphed
line_data <- Data %>% 
  group_by(Isolate, Readings,Dose_binary) %>% summarise(Mean_Abs = mean(Abs))
 #list of x intercepts e.g hour target is hit    
Xgrowth_Intercept<-line_data%>% group_by(Isolate,Dose_binary)%>% summarise(xint=(find_intersect(x=Readings,y=Mean_Abs,target = 0.2))/2)




# Multiple line plot
ggplot(Data, aes(
  x = (Readings / 2),
  y = Abs,
  color = Isolate
)) +
  geom_point(aes(color = Isolate)) +
  geom_vline(aes(xintercept = 24),
             linetype = "dashed",
             linewidth = 1.5) +
  geom_vline(aes(xintercept=xint),data = Xgrowth_Intercept)+
  scale_x_continuous(breaks = seq(0, 48, by = 12)) +
  scale_y_continuous(limits = c(0,1))+
  labs(y = "Percent Reduction of PrestoBlue", x = "Time (hours)") +
  facet_grid(~ Isolate + Dose_binary, scales = "free_y") +
  
  # coord_cartesian(xlim = c(0, 96))+
  #stat_smooth(method="loess", se=T, span = .5)+
  #scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal()




#math time!!!







#testing barbar

#subset to just one data 


DataA1 <- subset(Data,(Well %in% "A1"))

DataA1<-DataA1[,c(2,4)]


DataA2 <- subset(Data,(Well %in% "A2"))

DataA2<-DataA2[,c(2,4)]

DataA3 <- subset(Data,(Well %in% "A3"))

DataA3<-DataA3[,c(2,4)]




results_Bar4par <- Bayesfit(DataA1,model="Bar6par",inf.sigma=TRUE)

plot(DataA1,ylim=c(0,1))
t <- results_Bar4par$fit.t

y <- results_Bar4par$fit.ymean
lines(t,y,col="black",lwd=2)


t
y
summary(results_Bar4par)



set.seed(11)
results_H1 <- Bayescompare(DataA1, DataA2, hyp = "H1", model = "Bar4par")
e


log(10)*.2
