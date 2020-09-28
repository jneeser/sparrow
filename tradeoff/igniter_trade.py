import tradeoff_class as tc

# Categories
rest = tc.param(name="restartability",weight=0.084,Limitype ="fixed",Limit_val=[0,1])
comp = tc.param(name="system complexity",weight=0.224,Limitype ="fixed",Limit_val=[1,5])
mass = tc.param(name="system mass",weight=0.13,Limitype ="fixed",Limit_val=[1,5])
trl = tc.param(name="DARE TRL",weight=0.291,Limitype ="fixed",Limit_val=[1,9])
prob_ign = tc.param(name="probability of unintentional ignition",weight=0.116,Limitype ="fixed",Limit_val=[0,1])
imp_ign = tc.param(name="impact of unintentional ignition",weight=0.12,Limitype ="fixed",Limit_val=[1,5])
cost = tc.param(name="manufacturing cost",weight=0.035,Limitype ="fixed",Limit_val=[1,5])


# Designs 
pyro = tc.design(name="Pyrotechnic Igniter",sourcelist=[0,5,5,8,0,4,4.6666])
grain = tc.design(name="APCP igniter grain",sourcelist=[0,4,4.8333,5,0,2.8333,4])
ABS = tc.design(name="ABS Spark igniter",sourcelist=[1,2.333,3,3,1,4.16666,3.16666])
spark = tc.design(name="GOX/GH2 Sparktorch",sourcelist=[1,1.333,1.5,5,1,4.5,1.83333])

colors = [tc.color("EF5350", "red"), tc.color("FB8C00", "orange"), tc.color("FFEB3B", "yellow"), tc.color("8BC34A", "green"), tc.color("00BCD4", "blue")]

tradeoff_att =tc.tradeoff(design_list = [pyro,grain,ABS,spark],param_list= [rest,comp,mass,trl,prob_ign,imp_ign,cost])

tradeoff_att.get_tradeoff()
tradeoff_att.get_output(language="python",color_list=colors,width=10,rot="hor",caption="Tradeoff Ignition System")
#tradeoff_att.get_output()



sens = tc.sensitivity(tradeoff_att,20000)
sens.addto_technical(0.5)
sens.get_sens_linux()
print(sens.per)