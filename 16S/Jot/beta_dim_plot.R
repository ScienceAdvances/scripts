beta_dim_plot <- function(physeq,hue, method="PCoA",distance="bray",xlab,ylab,title,outdir,filename,w=7,h=6){
    ordination <- phyloseq::ordinate(physeq, method = method, distance = distance)
    px <- phyloseq::plot_ordination(physeq=physeq, ordination=ordination, type="samples", color = "Group")
    colnames(px$data)[1:2] <- c("V1","V2")
    data <- px$data
    if(method=="PCoA"){
        xlab <- paste(xlab,str_remove_all(px$labels$x,"(Axis.[12])|( )|(\\[)|(\\])"))
        ylab <- paste(ylab,str_remove_all(px$labels$y,"(Axis.[12])|( )|(\\[)|(\\])"))
    }
    p <- ggplot(data, aes(V1,V2)) +
        geom_point(size=8,aes(shape=Group,color=Group))+     #按分组画点的颜色和形状
        geom_text(aes(label=rownames(data)), size=4, vjust=-1, hjust=0)+  #添加样本名称,名称位置大小
        labs(x=xlab,y=ylab,title= title)+       #添加x、y轴标题和图标题
        geom_hline(yintercept=0, linetype=4) +        #添加原点垂直辅助线
        geom_vline(xintercept=0 ,linetype=4)+         #添加原点水平辅助线
        ggplot2::xlim(min(data$V1)*1.07,max(data$V1)*1.07)+
        ggplot2::ylim(min(data$V2)*1.07,max(data$V2)*1.07)+
        scale_color_manual(values=hue) +
        theme_bw(base_size = 18)+
        theme(plot.title=element_text(hjust=0.5))+setheme()
    using::gs(p,name=filename,outdir = outdir,w=w,h=h)

    p <- ggplot(data, aes(V1,V2)) +
        geom_point(size=8,aes(shape=Group,color=Group))+     #按分组画点的颜色和形状
        # geom_text(aes(label=rownames(data)), size=4, vjust=-1, hjust=0)+  #添加样本名称,名称位置大小
        labs(x=xlab,y=ylab,title= title)+       #添加x、y轴标题和图标题
        geom_hline(yintercept=0, linetype=4) +        #添加原点垂直辅助线
        geom_vline(xintercept=0 ,linetype=4)+         #添加原点水平辅助线
        ggplot2::xlim(min(data$V1)*1.07,max(data$V1)*1.07)+
        ggplot2::ylim(min(data$V2)*1.07,max(data$V2)*1.07)+
        scale_color_manual(values=hue) +
        theme_bw(base_size = 18)+
        theme(plot.title=element_text(hjust=0.5))+setheme()
    using::gs(p,name=paste0(filename,'_Lable'),outdir = outdir,w=w,h=h)

}