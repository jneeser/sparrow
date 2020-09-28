import tradeoff_class as tc

# Categories
cost = tc.param(name="Cost of external services",weight=0.155,Limitype ="fixed",Limit_val=[1,5])
isp = tc.param(name="Impact on ISP",weight=0.0863,Limitype ="SD",Limit_val=1.2)
mass = tc.param(name="system mass",weight=0.131,Limitype ="SD",Limit_val=1.4, direc="LB")
parts = tc.param(name="Fraction of reusable parts",weight=0.213,Limitype ="fixed",Limit_val=[1,5])
dev_comp = tc.param(name="Development complexity",weight=0.415,Limitype ="fixed",Limit_val=[1,5])



# Designs 
film = tc.design(name="Film Cooling",sourcelist=[4.333, 96, 10.56, 4.6667, 2.6667])
ablative = tc.design(name="Ablative Cooling",sourcelist=[4.333, 100, 6.27, 2, 4.3333])
regen = tc.design(name="Regenerative Cooling",sourcelist=[1.6667, 100, 0, 5, 2.333])
film_abl = tc.design(name="Ablative Cooling with Film Cooling",sourcelist=[3, 98, 8.42, 2.6667, 2.333])

colors = [tc.color("EF5350", "red"), tc.color("FB8C00", "orange"), tc.color("FFEB3B", "yellow"), tc.color("8BC34A", "green"), tc.color("00BCD4", "blue")]

tradeoff_att =tc.tradeoff(design_list = [film,ablative,regen,film_abl],param_list= [cost,isp,mass,parts,dev_comp])

tradeoff_att.get_tradeoff()
tradeoff_att.get_output(language="python",color_list=colors,width=10,rot="hor",caption="Tradeoff Cooling System")
#tradeoff_att.get_output()

sens = tc.sensitivity(tradeoff_att,200000)
sens.addto_technical(0.5)
sens.get_sens_linux()
print(sens.per)