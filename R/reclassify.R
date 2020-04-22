#' Make Powerpoint list for reclassification
#' @param fit1,fit2 An object of class glm
#' @param labels Optioanl labels for initial model and updated model
#' @param cutoff,cutoff2 A numeric vector Cutoff values for risk categories for initial and updated model
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' labels=c("Without Genetic Factor","With Genetic Factor")
#' data=makePPTList_reclassify(fit1,fit2,labels=labels)
makePPTList_reclassify=function(fit1,fit2,labels,cutoff=c(0,.1,.35,1),cutoff2){

   # cutoff=c(0,.1,.35,1)
  if(missing(labels)) labels=c("Initial Model","Updated Model")
  if(missing(cutoff2)) cutoff2=cutoff

  type=c("Rcode","Rcode","Rcode","table",rep("ggplot",5),"plot")
  title=c(labels[1],labels[2],
          "Reclassfication","Reclassfication Table",
          "Calibration Plot","Boxplot","predictPlot","priorPosteriorPlot",
          "riskDistributionPlot","ROC plot")
  temp=paste0("result<-reclassify(fit1,fit2, cutoff=",paste0("c(",paste0(cutoff,collapse=","),")"),
              ",cutoff2=",paste0("c(",paste0(cutoff,collapse=","),")"),");result")
  code=c(paste0("fit1<-",fit2call(fit1),";ORmultivariate(fit1)"),
         paste0("fit2<-",fit2call(fit2),";ORmultivariate(fit2)"),
         temp,
         paste0("reclassTable(result,labels=c('",paste0(labels,collapse="','"),"'))"),
         "calibrationPlot(fit2)","boxPlot2(fit2)",
         paste0("predictPlot(list(fit1,fit2),labels=c('",paste0(labels,collapse="','"),"'))"),
         "priorPosteriorPlot(list(fit1,fit2),alpha=0.3,color='red')",
         "riskDistributionPlot(fit2)",
         paste0("plotROC2(fit1,fit2,labels=c('",paste0(labels,collapse="','"),"'))"))
  data.frame(title=title,type=type,code=code,stringsAsFactors = FALSE)

}


#' Get call as string from an object of class glm
#' @param fit an obkect of class glm
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2call(fit1)
fit2call=function (fit) {
  temp = as.character(fit$call)
  tempcall=gsub(" ","",paste0(deparse(fit$formula),collapse = ""))
  footer = paste0(temp[1], "(", tempcall, ",family='", temp[3],
                  "',data=", temp[4])
  if (length(temp) > 4)
    footer = paste0(footer, ",offset=", temp[5])
  footer = paste0(footer, ")")
  footer
}



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
#'@examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' ORmultivariate(fit1)
#' ORmultivariate(fit2)
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
#' @param fit An object of class glm
#' @param ... further arguments to be passed to geom_point()
#' @importFrom Hmisc cut2
#' @importFrom ggplot2 geom_segment xlim ylim
#' @importFrom stats pchisq
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' calibrationPlot(fit2)
calibrationPlot=function(fit,...){

    y<-fit$y
    pred<-predict(fit,type="response")
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
#' @param fit An object of class glm
#' @param labels character optional labels
#' @param ... further arguments to be passed to geom_boxplot()
#' @importFrom ggplot2 geom_boxplot stat_summary
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' boxPlot2(fit2)
boxPlot2=function(fit,labels,...){
    y<-fit$y
    pred<-predict(fit,type="response")
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
#' @param fitlist A list of objects of class glm
#' @param labels String Optional labels
#' @param size numeric line size
#' @param ... Further arguments to be passed to geom_line()
#' @importFrom purrr map2_dfr map
#' @importFrom ggplot2 element_rect
#' @return A ggplot
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' labels=c("without genetic factors", "with genetic factors")
#' predictPlot(list(fit1,fit2),labels=labels)
predictPlot=function(fitlist,labels,size,...){
  if(missing(labels)) labels=c("Initial Model","Updated Model")
  predlist=map(fitlist,~predict(.x,type="response"))
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

    ggplot(data=df,aes_string(x="x",y="y",color="no",lty="no"))+geom_line(size=size,...)+
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
#' @param fitlist A list of objects of class glm
#' @param xlab,ylab,title character
#' @param ... Further arguments to be passed to geom_point()
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' priorPosteriorPlot(list(fit1,fit2),alpha=0.3,color="red")
priorPosteriorPlot=function(fitlist,xlab,ylab,title,...){
    if(missing(xlab)) xlab="Prior Risk"
    if(missing(ylab)) ylab="Posterior Risk"
    if(missing(title)) title="Prior versus Posterior Risk"

    predlist=map(fitlist,~predict(.x,type="response"))
    ggplot(data=NULL)+
        geom_point(aes(x=predlist[[1]],y=predlist[[2]]),...)+
        geom_segment(aes(x=0,y=0,xend=1,yend=1),color="black")+
        labs(x=xlab,y=ylab,title=title)+
        xlim(c(0,1))+ylim(c(0,1))+theme_bw()+
        theme_bw(base_size = 14)+
        theme(panel.grid = element_blank())

}



#' Risk Distribution Plot
#' @param fit An object of clas glm
#' @param labels Character
#' @param ... further arguments to be passed to geom_col()
#' @importFrom dplyr group_by mutate
#' @importFrom ggplot2 geom_col element_text element_rect
#' @importFrom rlang .data
#' @return A ggplot
#' @export
#' @examples
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' riskDistributionPlot(fit2)
riskDistributionPlot=function(fit,labels,...){
        # y=ExampleData$AMD;pred=pred2
        # labels=c("Without outcome", "With outcome")
    if(missing(labels)) labels=paste0(names(fit$model)[1]," = ", c(0,1))
    y=fit$y
    pred=predict(fit,type="response")
    m <- table(y, cut(pred, seq(0,1.05, 0.05)))
    df=as.data.frame(t(m))
    names(df)=c("x","no","Freq")
    df$x=as.numeric(gsub("\\(|,.*\\]","",df$x))

    df<-df %>% group_by(.data[["no"]]) %>% mutate(ratio=.data[["Freq"]]/sum(.data[["Freq"]]))
    df$no=factor(df$no,labels=labels)
    ggplot(df,aes_string(x="x",y="ratio",fill="no"))+
        geom_col(position="dodge")+
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
#' @param fit1,fit2 An object of class glm
#' @param labels Labels given to the ROC curves.
#' @importFrom PredictABEL plotROC
#' @importFrom pROC roc.test roc ci
#' @importFrom graphics text
#' @export
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' labels <- c("without genetic factors", "with genetic factors")
#' plotROC2(fit1,fit2,labels=labels)
plotROC2=function(fit1,fit2,labels){

    if(missing(labels)) labels <- c("Model 1", "Model 2")
    # data=ExampleData; cOutcome=2; predrisk=cbind(pred1,pred2); labels=labels
    data=fit1$data
    y=fit1$y
    cOutcome=which(names(fit1$data)==names(fit1$model)[1])
    pred1=predict(fit1,type="response")
    pred2=predict(fit2,type="response")
    predrisk=cbind(pred1,pred2)
    plotROC(data=data, cOutcome=cOutcome, predrisk=predrisk, labels=labels)
    roc1=roc(y,pred1)
    roc2=roc(y,pred2)
    auclabel=function(x){
        paste0("AUC: ",sprintf("%4.3f",x$auc),"(",sprintf("%4.3f",ci(x)[1]),"-",sprintf("%4.3f",ci(x)[3]),")") }
    auclabels=c(auclabel(roc1),auclabel(roc2))

    res=roc.test(roc1,roc2)
    delong="DeLong's test: "
    if(res$p.value<0.001) {
        delong=paste0(delong,"p < 0.001")
    } else {
        delong=paste0(delong,"p =", round(res$p.value,3))
    }
    text(0.75,0.25,label=delong)
    text(0.75,0.20,label=auclabel(roc1))
    text(0.75,0.15,label=auclabel(roc2),col="red")
}




#' Function for reclassification table and statistics
#' @param fit1,fit2 An object of class glm
#' @param cutoff Cutoff values for risk categories for initial model. Default value is c(0,.1,.35,1)
#' @param cutoff2 Optional cutoff values for risk categories for updated model
#' @importFrom Hmisc improveProb
#' @importFrom stats pnorm
#' @export
#' @return A list of class reclassified with elements
#' \describe{
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
#' }
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' result=reclassify(fit1,fit2)
#' result
reclassify=function (fit1,fit2, cutoff,cutoff2)
{

    pred1=predict(fit1,type="response")
    pred2=predict(fit2,type="response")
    y=fit1$y
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


#' Make reclassification table
#' @param x An object of class reclassified
#' @param outcome Character
#' @param labels Optional character
#' @importFrom officer fp_border
#' @importFrom flextable as_flextable align bold merge_v delete_part hline hline_top merge_at
#' add_footer_row as_grouped_data width add_header_row colformat_num hline_bottom
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @export
#' @return An object of class flextable
#' @examples
#' form1=paste0("AMD~",paste0(colnames(ExampleData)[3:10],collapse="+"))
#' form2=paste0("AMD~",paste0(colnames(ExampleData)[3:16],collapse="+"))
#' fit1=glm(as.formula(form1),data=ExampleData,family=binomial)
#' fit2=glm(as.formula(form2),data=ExampleData,family=binomial)
#' x=reclassify(fit1,fit2,cutoff=c(0,0.1,0.3,0.5,1))
#' reclassTable(x)
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
    a$inc=c(NA,inc,rep(NA,no-2))
    a$dec=c(NA,dec,rep(NA,no-2))
    a$space1=""
    a$nri=c(NA,x$x$nri.ne*100,rep(NA,no-2))
    a[[outcome]]="Absent"
    a
    b=x$present[,-ncol(x$present)]
    inc=sum(b[upper.tri(b)])
    dec=sum(b[lower.tri(b)])
    b=as.data.frame(b)
    b$space=""
    b$inc=c(NA,inc,rep(NA,no-2))
    b$dec=c(NA,dec,rep(NA,no-2))
    b$space1=""
    b$nri=c(NA,x$x$nri.e*100,rep(NA,no-2))
    b[[outcome]]="Present"
    df=rbind(a,b)
    x$cutoff2

    colnames(df)[1:(length(x$cutoff2)-1)]=cutoff2str(x$cutoff2)
    df[["initial"]]=rep(cutoff2str(x$cutoff),2)
    names(df)[names(df)=="space"]=" "
    names(df)[names(df)=="space1"]="1"
    #str(df)
    df<-df %>% select(all_of("initial"),everything())
    df<-as_grouped_data(df,groups=outcome)
    df

    ft<-as_flextable(df) %>% bold(j=1,i=~!is.na(outcome),bold=TRUE,part="body")


    ft<-delete_part(ft,part="header")
    j=c(cutoff2str(x$cutoff2),"inc","dec")
    ft<-colformat_num(ft,j=j,big.mark=",",digits=0,na_str="")
    ft<-colformat_num(ft,j="nri",big.mark=",",digits=1,na_str="")

    space=no+2
    space1=no+5
    ft<-width(ft,j=space,width=0.1)
    ft<-width(ft,j=space1,width=0.1)
    values=c("",cutoff2str(x$cutoff2),"","Increased risk","Decreased risk","","")
    ft<-add_header_row(ft,values=values)
    big_border = fp_border(color="black", width = 2)
    small_border = fp_border(color="black", width = 1)
    ft<-hline_bottom(ft,border=big_border,part="header")
    values=c(labels[1],labels[2],"","Reclassified","","Net correctly reclassified(%)")
    colwidths=c(1,no,1,2,1,1)
    ft<-add_header_row(ft,values=values,colwidths=colwidths)
    ft<-width(ft,j=2:(no+1),width=1.4)

    ft<-hline_top(ft,part="body",border=big_border)
    ft<-hline(ft,i=1,j=2:(no+1),part="header",border=small_border)
    ft<-hline(ft,i=1,j=(no+3):(no+4),part="header",border=small_border)
    ft<-hline_top(ft,part="header",border=big_border)
    ft
    ft<-merge_at(ft,i=1:2,j=space1+1,part="header")
    ft<-merge_at(ft,i=1:2,j=1,part="header")
    ft<-width(ft,j=(space+1):(space+2),width=1.2)
    ft<-width(ft,j=space1+1,width=0.6)

    y=x$NRIcat
    values=paste0("Net reclassification improvement :",round(y[1]*100,1),"% (95% CI :",
                  round(y[2]*100, 1), "-", round(y[3]*100, 1),"); p ")

    if(y[4]<0.001) {
        values=paste0(values,"< .001")
    } else {
        values=paste0(values,"=", round(y[4],3))
    }
    ft<-add_footer_row(ft,values=values,colwidths=space1+1)
    ft<-align(ft,align="center",part="header")
    y=x$IDI
    values=paste0("Integrated discrimination improvement :",round(y[1]*100,1),"% (95% CI :",
                  round(y[2]*100, 1), "-", round(y[3]*100, 1),"); p ")
    if(y[4]<0.001) {
        values=paste0(values,"< .001")
    } else {
        values=paste0(values,"=", round(y[4],3))
    }
    ft<-add_footer_row(ft,top=FALSE,values=values,colwidths=space1+1)
    ft<-flextable::fontsize(ft,size=12,part="all")
    ft<-hline(ft,i=1,j=1,part="header",border=big_border)
    ft<-hline(ft,i=1,j=space1+1,part="header",border=big_border)
    ft
}


#' Convert cutoff to string
#' @param cutoff Numeric vector
#' @export
#' @examples
#' cutoff2str(cutoff=c(0,.1,.2,.5,1))
cutoff2str=function(cutoff=c(0,.1,.35,1)){
    # cutoff=c(0,.1,.35,1)
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

