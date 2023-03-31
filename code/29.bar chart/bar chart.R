library(tidyverse)
library(rstatix)    # 整洁统计检验

df <- test
head(df)
df1 <- df%>%
  pivot_longer(-1,names_to="lncRNA",values_to="value")
df1

df1%>%
  group_by(risk,lncRNA)%>%
  shapiro.test(value)

ggplot(df1,aes(risk=value))+
  stat_qq()+
  stat_qq_line(color="red")+
  facet_wrap(risk~lncRNA,scales="free")

df2=df1%>%
  group_by(lncRNA,risk)%>%
  summarise(Relative_expression_of_lncRNA=mean(value),se=sd(value))
df2

df2%>%
  ggplot(aes(lncRNA,Relative_expression_of_lncRNA,fill=risk))+
  geom_col(position="dodge")+
  geom_errorbar(aes(ymax=Relative_expression_of_lncRNA+se,ymin=Relative_expression_of_lncRNA-se),width=0.15,
                position=position_dodge(0.9))+
  ylim(0,3.7)+
  annotate("text",x=1,y=2+1.5,label="**")+#
  annotate("text",x=2,y=2+1.5,label="***")+
  annotate("text",x=3,y=2+1.5,label="*")+
  annotate("text",x=4,y=2+1.5,label="**")+
  annotate("text",x=5,y=2+1.5,label="**")+
  annotate("text",x=6,y=2+1.5,label="*")+
  annotate("text",x=7,y=2+1.5,label="*")+
  annotate("text",x=8,y=2+1.5,label="***")

