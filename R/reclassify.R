# library(PredictABEL)
#
# ### reclassification
# data(ExampleData)
# names(ExampleData)[2]="AMD"
# usethis::use_data(ExampleData)
#
# form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
# form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
# fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
# fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
# pred1=predict(fit1,type="response")
# pred2=predict(fit2,type="response")
#
#
#
# require(ggplot2)
# require(purrr)
# require(dplyr)
# require(pROC)
# calibrationPlot(y=ExampleData$AMD,pred2)
# boxPlot2(y=ExampleData$AMD,pred2)
# labels=c("without genetic factors", "with genetic factors")
# predictPlot(list(pred1,pred2),labels=labels)
# priorPosteriorPlot(pred1,pred2,alpha=0.3,color="red")
# riskDistributionPlot(y=ExampleData$AMD,pred=pred2)
# labels <- c("without genetic factors", "with genetic factors")
# plotROC2(data=ExampleData, cOutcome=2, predrisk=cbind(pred1,pred2), labels=labels)
# plotROC2(data=ExampleData, cOutcome=2, predrisk=cbind(pred1,pred2))
# reclassification(data=ExampleData, cOutcome=2,
#                  predrisk1=pred1, predrisk2=pred2, cutoff <- c(0,.1,.35,1))


#'A hypothetical dataset that is used to demonstrate all functions.
#'ExampleData is a hypothetical dataset constructed to demonstrate all functions in the package. ExampleData is a data frame containing a binary outcome variable (e.g., disease present/absent) and genetic and non-genetic predictor variables for 10,000 persons. In this dataset, column 1 is the ID number, column 2 is the outcome variable, columns 3-10 are non-genetic variables and columns 11-16 are genetic variables.
#'@examples
#'data(ExampleData)
#'head(ExampleData,5)
#'@references
#'Suman Kundu, Yurii S. Aulchenko and A. Cecile J.W. Janssens (2020). PredictABEL: Assessment of Risk Prediction Models. R package version 1.2-4.https://CRAN.R-project.org/package=PredictABEL
"ExampleData"


#'Function to obtain multivariate odds ratios from a logistic regression model.
#'@param fit An object of class glm
#'@importFrom methods is
#'@export
#'@return A list with elements
#'\describe{
#' \item{Predictors_Summary}{OR with 95\% CI and corresponding p-values for each predictor in the model}
#' \item{Brier_Score}{Brier score}
#' \item{Nagelkerke_Index}{Nagelkerke's R2 value}
#'}
#'@references
#'Suman Kundu, Yurii S. Aulchenko and A. Cecile J.W. Janssens (2020). PredictABEL: Assessment of Risk Prediction Models. R package version 1.2-4.https://CRAN.R-project.org/package=PredictABEL
#'Brier GW. Verification of forecasts expressed in terms of probability. Monthly weather review 1950;78:1-3.
#'Nagelkerke NJ. A note on a general definition of the coefficient of determination. Biometrika 1991;78:691-692.
ORmultivariate=function (fit)
{
  if (!is(fit, "glm"))
    stop("fit argument should have 'glm' class")
  sum.coef <- summary(fit)$coef
  OR <- exp(sum.coef[, "Estimate"])
  Upper.CI <- exp(sum.coef[, "Estimate"] + 1.96 * sum.coef[,
                                                           "Std. Error"])
  Lower.CI <- exp(sum.coef[, 1] - 1.96 * sum.coef[, 2])
  tab <- cbind(OR = round(OR, 4), Lower.CI = round(Lower.CI,
                                                   4), Upper.CI = round(Upper.CI, 4), `p-value` = round((sum.coef)[,
                                                                                                                   4], 4))
  B <- mean((fit$y) * (1 - predict(fit, type = "response"))^2 +
              (1 - fit$y) * (predict(fit, type = "response"))^2)
  Bmax <- mean(fit$y) * (1 - mean(fit$y))^2 + (1 -
                                                             mean(fit$y)) * mean(fit$y)^2
  Bscaled <- 1 - B/Bmax
  LLfull <- (fit$deviance)/(-2)
  LLnul <- (fit$null.deviance)/(-2)
  n <- length(fit$y)
  Mc <- 1 - (LLfull/LLnul)
  Cox <- 1 - exp((-2) * (LLfull - LLnul)/n)
  Rmax <- 1 - exp((2 * LLnul)/n)
  Nag <- Cox/Rmax
  p <- list(Predictors_Summary = tab, Brier_Score = round(B,4),
            Nagelkerke_R2 = round(Nag, 4))
  p
}


#' Draw a calibration plot
#' @param y A numeric vector
#' @param pred Numeric vecor as a result of predict.glm
#' @param ... further arguments to be passed to geom_point()
#' @importFrom Hmisc cut2
#' @importFrom ggplot2 geom_segment xlim ylim
#' @importFrom stats pchisq
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred2=predict(fit2,type="response")
#' calibrationPlot(y=ExampleData$AMD,pred2)
calibrationPlot=function(y,pred,...){

    groups=10
    p<-pred
    matres <- matrix(NA, nrow = groups, ncol = 5)
    sor <- order(p)
    p <- p[sor]
    y <- y[sor]
    groep <- Hmisc::cut2(p, g = groups)
    total <- tapply(y, groep, length)
    predicted <- round(tapply(p, groep, sum), 2)
    observed <- tapply(y, groep, sum)
    meanpred <- round(tapply(p, groep, mean), 3)
    meanobs <- round(tapply(y, groep, mean), 3)
    matres <- cbind(total, meanpred, meanobs, predicted,
                    observed)
    contr <- ((observed - predicted)^2)/(total * meanpred *
                                             (1 - meanpred))
    chisqr <- sum(contr)
    df <- (groups - 2)
    pval <- 1 - pchisq(chisqr, df)
    temp=paste0("atop(chi^2 == ",round(chisqr,2),"(df == ",df,"),italic(p) ")
    if(pval<0.001){
        temp=paste0(temp,"< 0.001)")
    }else {
        temp=paste0(temp,"== ", pval,")")
    }

    df1=as.data.frame(matres)

    ggplot(data=df1,aes_string(x="meanpred",y="meanobs"))+
        geom_segment(x=0,y=0,xend=1,yend=1,color="grey50")+
        geom_point(...)+
        labs(x="Predicted risk",y="Observed risk",title="Calibration plot")+
        xlim(c(0,1))+ylim(c(0,1))+
        annotate("text",x=0.65,y=0.2,label=temp,hjust=0,parse=TRUE)+
        theme_bw()+
        theme(panel.grid = element_blank())
}


#' Draw a box plot
#' @param y A numeric vector
#' @param pred Numeric vecor as a result of predict.glm
#' @param ... further arguments to be passed to geom_boxplot()
#' @importFrom ggplot2 geom_boxplot stat_summary
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred2=predict(fit2,type="response")
#' boxPlot2(y=ExampleData$AMD,pred=pred2)
boxPlot2=function(y,pred,labels,...){
    if(missing(labels)) labels=c("Without disease", "With disease")
    df=data.frame(y=pred,x=y)
    meandiff=mean(df$y[df$x==unique(df$x)[2]])-mean(df$y[df$x==unique(df$x)[1]])
    df$x=factor(df$x,labels=labels)
    ggplot(data=df,aes_string(x="x",y="y",color="x"))+
        geom_boxplot(...)+
        stat_summary(fun="mean",geom="point",pch=15,size=4)+
        labs(x="",y="Predicted risks",title="Difference between groups")+
        annotate("text",x=1.5,y=0.8,
                 label=paste0("Mean diff = ",round(meandiff,3)))+
        theme_bw(base_size = 14)+
        theme(panel.grid = element_blank(),legend.position="None")
}


#' Draw a Predictiveness Plot
#' @param predlist A list of prediction vectors
#' @param labels String Optional labels
#' @param ... Further arguments to be passed to geom_line()
#' @return A ggplot
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' labels=c("without genetic factors", "with genetic factors")
#' predictPlot(list(pred1,pred2),labels=labels)
predictPlot=function(predlist,labels,size,...){
    count=length(predlist)
    if(missing(labels)) labels <- paste0("Model",1:count)
    if(missing(size)) size=1
    df<-map2_dfr(predlist,1:count,function(pred,i){
        y <- sort(pred)
        n <- length(y)
        x <- (1:n)/n
        df=data.frame(x=x,y=y)
        df$no=i
        df
    })


    df$no=factor(df$no,labels=labels)

    ggplot(data=df,aes(x=x,y=y,color=no,lty=no))+geom_line(size=size,...)+
        labs(x="Cumulative percentage",
             y="Predicted risks",
             title="Predictiveness Curve")+
        theme_bw(base_size = 14)+
        theme(legend.position=c(0.3,0.9),
              legend.title = element_blank(),
              legend.box.background = element_rect(),
              panel.grid = element_blank())
}


#' Draw prior versus posterior plot
#' @param pred1 Numeric vecor as a result of predict.glm
#' @param pred2 Numeric vecor as a result of predict.glm
#' @param xlab,yla,title character
#' @param ... Further arguments to be passed to geom_point()
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' priorPosteriorPlot(pred1,pred2,alpha=0.3,color="red")
priorPosteriorPlot=function(pred1,pred2,xlab,ylab,title,...){
    if(missing(xlab)) xlab="Prior Risk"
    if(missing(ylab)) ylab="Posterior Risk"
    if(missing(title)) title="Prior versus Posterior Risk"
    ggplot(data=NULL)+
        geom_segment(aes(x=0,y=0,xend=1,yend=1),color="black")+
        geom_point(aes(x=pred1,y=pred2),...)+
        labs(x=xlab,y=ylab,title=title)+
        xlim(c(0,1))+ylim(c(0,1))+theme_bw()+
        theme_bw(base_size = 14)+
        theme(panel.grid = element_blank())

}



#' Risk Distribution Plot
#' @param y A numeric vector
#' @param pred Numeric vecor as a result of predict.glm
#' @param labels Character
#' @param ... further arguments to be passed to geom_col()
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred2=predict(fit2,type="response")
#' riskDistributionPlot(y=ExampleData$AMD,pred=pred2)
riskDistributionPlot=function(y,pred,labels,...){
    if(missing(labels)) labels=c("Without outcome", "With outcome")

    m <- table(y, cut(pred, seq(0,1.05, 0.05)))
    df=as.data.frame(t(m))
    names(df)=c("x","no","Freq")
    df$x=as.numeric(gsub("\\(|,.*\\]","",df$x))
    df<-df %>% group_by(no) %>% mutate(ratio=Freq/sum(Freq))
    df$no=factor(df$no,labels=labels)
    ggplot(df,aes(x=x,y=ratio,fill=no))+
        geom_col(position="dodge",...)+
        labs(x="Predicted risk",y="Precentage",
             title="Distribution of Predicted Risk")+
        theme_bw()+
        theme(legend.position=c(0.8,0.9),
              legend.title = element_blank(),
              legend.box.background = element_rect(),
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5))
}





#' Draw ROC curves
#' @param data	A data.frame
#' @param cOutcome Column number of the outcome variable.
#' @param predrisk	Vector of predicted risk
#' @param labels Labels given to the ROC curves.
#' @importFrom PredictABEL plotROC
#' @importFrom pROC roc.test
#' @export
#' @examples
#' labels <- c("without genetic factors", "with genetic factors")
#' plotROC2(data=ExampleData, cOutcome=2, predrisk=cbind(pred1,pred2), labels=labels)
plotROC2=function(data, cOutcome, predrisk, labels){

    if(missing(labels)) labels <- c("Model 1", "Model 2")
    plotROC(data=data, cOutcome=cOutcome, predrisk=predrisk, labels=labels)
    roc1=roc(data[[cOutcome]],predrisk[1])
    roc2=roc(data[[cOutcome]],predrisk[2])
    auclabel=function(x){
        paste0("AUC: ",round(x$auc,3),"(",round(ci(x)[1],3),"-",round(ci(x)[3],3),")") }
    auclabels=c(auclabel(roc1),auclabel(roc2))

    res=roc.test(roc1,roc2)
    delong="DeLong's test: "
    if(res$p.value<0.001) {
        delong=paste0(delong,"p < 0.001")
    } else {
        delong=paste0(delong,"p =", round(res$p.value,3))
    }
    text(0.7,0.25,label=delong)
    text(0.7,0.20,label=auclabel(roc1))
    text(0.7,0.15,label=auclabel(roc2),col="red")
}




#' Draw ROC curves
#' @param y A numeric vector
#' @param predlist A list of numeric vecors as results of predict.glm
#' @param labels Character
#' @param ... further arguments to be passed to geom_line()
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' labels <- c("without genetic factors", "with genetic factors")
#' ggplotROC(ExampleData$AMD,predlist=list(pred1,pred2),labels=labels)
ggplotROC=function(y,predlist,labels,...) {

if(missing(labels)) labels <- c("Model 1", "Model 2")

count=length(predlist)
temp=list()
df<-predlist %>% map2_dfr(1:count,function(a,b){
    temp[[b]]<<-roc(y~a)
    df=data.frame(y=temp[[b]]$sensitivities,x=1-temp[[b]]$specificities)
    df$no=b
    df
})

auclabel=temp %>% map_chr(function(x){
    paste0("AUC: ",round(x$auc,3),"(",round(ci(x)[1],3),"-",round(ci(x)[3],3),")")
})
df2=data.frame(label=auclabel)
df2$no=1:count
df2$x=0.8
df2$y=0.3-(df2$no/20)
df2$no=factor(df2$no,labels=labels)
res=roc.test(temp[[1]],temp[[2]])

df$no=factor(df$no,labels=labels)
delong="DeLong's test: "
if(res$p.value<0.001) {
    delong=paste0(delong,"p < 0.001")
} else {
    delong=paste0(delong,"p =", round(res$p.value,3))
}
p<-ggplot(data=df)+
    geom_line(aes(x=x,y=y,color=no,group=no,lty=no),size=1,...)+
    labs(x="1-Specificity",y="Sensitivity",
         title="ROC plot")+
    geom_segment(x=0,y=0,xend=1,yend=1,color="grey50",lty=2)+
    theme_bw()+
    theme(legend.position=c(0.8,0.1),
          legend.title = element_blank(),
          legend.box.background = element_rect(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
p+geom_text(data=df2,aes(x=x,y=y,label=label,color=no))+
    annotate("text",x=0.8,y=0.3,label=delong)
}





#' Function for reclassification table and statistics
#' @param y A numeric vector of outcome
#' @param pred1,pred2 A numeric vector as a result of predict.glm
#' @param cutoff Cutoff values for risk categories for initial model. Default value is c(0,.1,.35,1)
#' @param cutoff2 Optional cutoff values for risk categories for updated model
#' @importFrom Hmisc improveProb
#' @export
#' @return A list of class reclassified with elements
#' \item{cutoff}{Cutoff values for risk categories for initial model}
#' \item{cutoff2}{Cutoff values for risk categories for updated model}
#' \item{absent}{summary table for absent outcome}
#' \item{present}{summary table for present outcome}
#' \item{combined}{summary table for combined outcome}
#' \item{x}{improveProb for categorical variables}
#' \item{y}{improveProb for initial and updated model}
#' \item{NRIcat}{Numeric vector for Net reclassification improvement, categorical}
#' \item{NRIcont}{Numeric vector for Net reclassification improvement, continuous}
#' \item{IDI}{Numeric vector for Integrated discrimination improvement}
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' result=reclassify(ExampleData$AMD,pred1,pred2)
#' result
reclassify=function (y, pred1, pred2, cutoff,cutoff2)
{
    if(missing(cutoff)) cutoff=c(0,.1,.35,1)
    if(missing(cutoff2)) cutoff2=cutoff

    c1 <- cut(pred1, breaks = cutoff, include.lowest = TRUE,
              right = FALSE)
    c2 <- cut(pred2, breaks = cutoff2, include.lowest = TRUE,
              right = FALSE)
    tabReclas <- table(`Initial Model` = c1, `Updated Model` = c2)
    ta <- table(c1, c2, y)
    TabAbs <- ta[, , 1]
    tab1 <- cbind(TabAbs, ` % reclassified` = round((rowSums(TabAbs) -
                                                         diag(TabAbs))/rowSums(TabAbs), 2) * 100)
    names(dimnames(tab1)) <- c("Initial Model", "Updated Model")

    TabPre <- ta[, , 2]
    tab2 <- cbind(TabPre, ` % reclassified` = round((rowSums(TabPre) -
                                                         diag(TabPre))/rowSums(TabPre), 2) * 100)
    names(dimnames(tab2)) <- c("Initial Model", "Updated Model")

    Tab <- tabReclas
    tab <- cbind(Tab, ` % reclassified` = round((rowSums(Tab) -
                                                     diag(Tab))/rowSums(Tab), 2) * 100)
    names(dimnames(tab)) <- c("Initial Model", "Updated Model")


    c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
    c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
    x <- Hmisc::improveProb(x1 = as.numeric(c11) * (1/(length(levels(c11)))),
                     x2 = as.numeric(c22) * (1/(length(levels(c22)))), y = y)

    y <- Hmisc::improveProb(x1 = pred1, x2 = pred2, y = y)

    NRIcat=c(x$nri,x$nri - 1.96 * x$se.nri,x$nri +1.96 * x$se.nri,2 * pnorm(-abs(x$z.nri)))

    NRIcont=c(y$nri,y$nri - 1.96 * y$se.nri,y$nri +1.96 * y$se.nri,2 * pnorm(-abs(y$z.nri)))

    IDI=c(y$idi,y$idi -1.96 * y$se.idi,y$idi + 1.96 * y$se.idi,2 * pnorm(-abs(y$z.idi)))


    result=list(
        cutoff=cutoff,
        cutoff2=cutoff2,
        absent=tab1,
        present=tab2,
        combined=tab,
        x=x,
        y=y,
        NRIcat=NRIcat,
        NRIcont=NRIcont,
        IDI=IDI
    )

    class(result)="reclassified"
    invisible(result)
}

#' S3 method print for an object of class reclassified
#' @param x An object of class reclassified
#' @param ... Fruther arguments
print.reclassified=function(x,...){
    cat(" _________________________________________\n")
    cat(" \n     Reclassification table    \n")
    cat(" _________________________________________\n")
    cat("\n Outcome: absent \n  \n")
    print(x$absent)
    cat("\n \n Outcome: present \n  \n")
    print(x$present)
    cat("\n \n Combined Data \n  \n")
    print(x$combined)
    cat(" _________________________________________\n")
    y=x$NRIcat
    cat("\n NRI(Categorical) [95% CI]:", round(y[1], 4), "[",
        round(y[2], 4), "-", round(y[3], 4), "]", "; p-value:", round(y[4], 5), "\n")
    y=x$NRIcont
    cat(" NRI(Continuous) [95% CI]:", round(y[1], 4), "[",
        round(y[2], 4), "-", round(y[3], 4), "]", "; p-value:", round(y[4], 5), "\n")
    y=x$IDI
    cat(" IDI [95% CI]:", round(y[1], 4), "[",
        round(y[2], 4), "-", round(y[3], 4), "]", "; p-value:", round(y[4], 5), "\n")
}


#' Make recalssification table
#' @param x An object of class reclassified
#' @param outcome Character
#' @param labels Optional character
#' @importFrom officer fp_border
#' @importFrom flextable as_flextable
#' @export
#' @return An object of class flextable
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' pred1=predict(fit1,type="response")
#' pred2=predict(fit2,type="response")
#' result=reclassify(ExampleData$AMD,pred1,pred2)
#' reclassTable(result)
reclassTable=function(x,outcome,labels){
    if(missing(outcome)) outcome="outcome"
    if(missing(labels)) labels=c("Initial Model","Updated Model")

    # outcome="outcome"
    # labels=c("Initial Model","Updated Model")

    a=x$absent[,-ncol(x$absent)]
    no=ncol(a)
    inc=sum(a[upper.tri(a)])
    dec=sum(a[lower.tri(a)])
    a=as.data.frame(a)
    a$space=""
    a$inc=inc
    a$dec=dec
    a$space1=""
    a$nri=x$x$nri.ne*100
    a[[outcome]]="Absent"
    a
    b=x$present[,-ncol(x$present)]
    inc=sum(b[upper.tri(b)])
    dec=sum(b[lower.tri(b)])
    b=as.data.frame(b)
    b$space=""
    b$inc=inc
    b$dec=dec
    b$space1=""
    b$nri=x$x$nri.e*100
    b[[outcome]]="Present"
    df=rbind(a,b)

    colnames(df)[1:(length(x$cutoff2)-1)]=cutoff2str(x$cutoff2)
    df$initial=rep(cutoff2str(x$cutoff),2)
    names(df)[names(df)=="space"]=" "
    names(df)[names(df)=="space1"]="1"
    str(df)
    df<-df %>% select(initial,everything())
    df<-as_grouped_data(df,groups=outcome)
    df
    ft<-as_flextable(df) %>% bold(j=1,i=~!is.na(outcome),bold=TRUE,part="body")

    ft<-merge_v(ft,j=c("inc","dec","nri"),part="body")

    ft<-delete_part(ft,part="header")
    j=c(cutoff2str(x$cutoff2),"inc","dec")
    ft<-colformat_num(ft,j=j,big.mark=",",digits=0)
    ft<-colformat_num(ft,j="nri",big.mark=",",digits=1)

    ft<-width(ft,j=5,width=0.1)
    ft<-width(ft,j=8,width=0.1)
    values=c("",cutoff2str(x$cutoff2),"","Increased risk","Decreased risk","","")
    ft<-add_header_row(ft,values=values)
    ft
    values=c(labels[1],labels[2],"","Reclassified","","Net correctly reclassified(%)")
    colwidths=c(1,no,1,2,1,1)
    ft<-add_header_row(ft,values=values,colwidths=colwidths)
    ft<-width(ft,j=2:4,width=1.2)
    big_border = fp_border(color="black", width = 2)
    small_border = fp_border(color="black", width = 1)
    ft<-hline_top(ft,part="body",border=big_border)
    ft<-hline(ft,i=1,j=2:4,part="header",border=small_border)
    ft<-hline(ft,i=1,j=6:7,part="header",border=small_border)
    ft<-hline_top(ft,part="header",border=big_border)
    ft
    ft<-merge_at(ft,i=1:2,j=9,part="header")
    ft<-merge_at(ft,i=1:2,j=1,part="header")
    ft<-width(ft,j=6:7,width=1.1)
    ft<-width(ft,j=9,width=0.5)

    y=x$NRIcat
    values=paste0("Net reclassification improvement :",round(y[1]*100,1),"% (95% CI :",
                  round(y[2]*100, 1), "-", round(y[3]*100, 1),"); p ")

    if(y[4]<0.001) {
        values=paste0(values,"< .001")
    } else {
        values=paste0(values,"=", round(y[4],3))
    }
    ft<-add_footer_row(ft,values=values,colwidths=9)
    ft<-align(ft,align="center",part="header")
    y=x$IDI
    values=paste0("Integrated discrimination improvement :",round(y[1]*100,1),"% (95% CI :",
                  round(y[2]*100, 1), "-", round(y[3]*100, 1),"); p ")
    if(y[4]<0.001) {
        values=paste0(values,"< .001")
    } else {
        values=paste0(values,"=", round(y[4],3))
    }
    ft<-add_footer_row(ft,top=FALSE,values=values,colwidths=9)
    ft
}


#' Convert cutoff to string
#' @param cutoff Numeric vector
#' @export
#' @examples
#' cutoff2str(cutoff=c(0,.1,.35,1))
cutoff2str=function(cutoff=c(0,.1,.35,1)){
    cutoff=c(0,.1,.35,1)
    no=length(cutoff)
    cutoff=cutoff*100
    temp=paste0("< ",cutoff[2])
    if(no>3){
    for(i in 2:(no-2)){
        temp=c(temp,paste0(cutoff[i],"-",cutoff[i+1]))
    }
    }
    temp=c(temp,paste0("\U2265 ",cutoff[no-1]))
    paste0(temp,"%")
}

