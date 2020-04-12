#' Draw histogram of transformed value with normal density
#' @param x A numeric vector
#' @param best logical Whether or not use bestNormalize::bestNormalize
#' @param power logical Whether or not use car::powerTransform
#' @param hist logical Whether or not draw histogram
#' @param density logical Whether or not draw density
#' @param select numeric choice of normality test
#' @param method character choice of normality test
#' @importFrom car powerTransform
#' @importFrom bestNormalize bestNormalize
#' @importFrom stats shapiro.test na.omit
#' @importFrom nortest ad.test cvm.test lillie.test pearson.test sf.test
#' @importFrom purrr map_df map2_df map_dfr `%>%`
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom ggplot2 ggplot geom_histogram geom_density facet_wrap geom_text theme theme_bw
#' @importFrom ggplot2 geom_line scale_x_continuous element_blank labs aes_string guide_axis
#' @export
#' @examples
#' require(moonBook)
#' transformDensity(acs$TG)
#' transformDensity(acs$TG,best=TRUE)
transformDensity=function(x,
                          best=FALSE,
                          power=TRUE,
                          hist=TRUE,
                          density=FALSE,
                          select=NULL,
                          method=c("shapiro","ad","cvm","lillie","pearson","sf")){

        xname=deparse(substitute(x))

        funs=c("shapiro","ad","cvm","lillie","pearson","sf")
        if(!missing(select)) {
                method=funs[select]
        } else{
                method=match.arg(method)
        }
        res=summary(powerTransform(x))
        rhambda=res$result[1,1]

        x=na.omit(x)
        if(best){
                x=bestNormalize(x,r = 1, k = 5, warn = F)
                result=list()
                colNames=list()
                result[[1]]=x$chosen_transform$x.t
                colNames[[1]]=attr(x$chosen_transform,"class")[1]

                for(i in 1:length(x$other_transforms)){
                        result[[i+1]]=x$other_transforms[[i]]$x.t
                        colNames[[i+1]]=attr(x$other_transforms[[i]],"class")[1]
                }
                df=as.data.frame(result)
                names(df)=unlist(colNames)
                df
        } else{
           levels=c("1/x^2","1/x","1/sqrt(x)","log(x)","sqrt(x)","x","x^2")
           if(power) levels=c(levels,paste0("x^",round(rhambda,3)))
           names(levels)=levels
           df<-map_df(levels,~eval(parse(text=.x)))
        }
        df2<-df %>% map_df(~seq(min(.x),max(.x),length.out=50))
        df3<-map2_df(df2,df,~dnorm(.x,mean=mean(.y),sd=sd(.y)))
        longdf2=pivot_longer(df2,cols=everything(),values_to="x")
        longdf3=pivot_longer(df3,cols=everything(),values_to="y")
        longdf4=cbind(longdf2,longdf3[,2])

        longdf=pivot_longer(df,cols=everything())
        longdf

        if(method=="shapiro"){
                fun=shapiro.test
        } else if(method=="ad") {
                fun=ad.test
        } else if(method=="cvm"){
                fun=cvm.test
        } else if(method=="lillie"){
                fun=lillie.test
        } else if(method=="pearson"){
                fun=pearson.test
        } else{
                fun=sf.test
        }

        p2character=function (x, digits = 3)
        {
                cut = 1/(10^digits)
                temp = sprintf(paste0("%.", digits, "f"), x)
                temp = gsub("^0", "",temp)

                for (i in 1:length(x)) {
                        if (is.na(x[i])) {
                                temp[i] = ""
                        }
                        else if (x[i] < cut) {
                                temp[i] = gsub("0$", "1",temp[i])
                                temp[i] = paste0("italic(p) < ", temp[i])
                        } else{
                                temp[i] = paste0("italic(p) == ", temp[i])
                        }
                }
                temp
        }

        if(best){
                df5=data.frame(name=names(x$norm_stats),value=round(x$norm_stats,3),stringsAsFactors = FALSE)
        } else{
                df5<-df %>% map_dfr(~fun(.x) %>%
                                            .$p.value %>%
                                            p2character) %>%
                        pivot_longer(cols=everything())
        }

        subtitle=ifelse(best,
                        "Estimated Normality Statistics (Pearson P/df: the lower, the more normal)",
                        paste0("p values by ",method,".test"))

        p<-ggplot(data=longdf,aes_string(x="value"))
        if(hist) p<-p+
                geom_histogram(aes_string(y="..density.."),color="grey50",fill="cornsilk",alpha=0.5)
        if(density) p<-p+geom_density()
        p<-p+   facet_wrap(~name,scales="free")+
                geom_text(data=df5,x=Inf,y=Inf,aes_string(label="value"),hjust=1.2,vjust=2,parse=TRUE)+
                theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
                geom_line(data=longdf4,aes_string(x="x",y="y"),color = "red")+
                scale_x_continuous(guide = guide_axis(n.dodge = 2))+
                labs(x="",title=paste0("Transformation of ",xname),subtitle=subtitle)
        p

}

