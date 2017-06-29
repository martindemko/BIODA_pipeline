#!/usr/bin/env python

# required module(s)
import json
# default and configuration files (may required full path)
default_configuration_file = "defaults.txt"
user_json_configuration_file = "config_2D_test.json"

#load user specific configuration json file
def load_user_config_file(json_conf_file): #json_conf_file is file with full path 
	with open(json_conf_file) as json_data: config = json.load(json_data)
	print "[$] user config json file loaded"
	return(config)

#define full path to default config file for snakemake:
#create dict of defaults config file for snakemake:
def load_default_config_file(def_conf_file):
	defaults={}; main_category="Na"; sub_category="Na";
	f=open(def_conf_file, "r")
	for line in f.xreadlines():
		if line[0]=="#": pass
		else:
			if line in ['\n', '\r\n']: main_category="Na"; sub_category="Na"; #detect new category
			else:
				try:
					column=line.strip().split("="); feature=column[0]; value=column[1]
					if feature == "category": main_category=value; defaults[main_category]={}
					else:
						if feature==main_category:
							sub_category=value; defaults[main_category][sub_category]={}
						else:
							if main_category=="Na": defaults[feature]=value
							elif sub_category != "Na": defaults[main_category][sub_category][feature]=value
							else: defaults[main_category][feature]=value
				except: pass
	f.close()
	print "[$] default config file loaded"
	return(defaults)

def update_user_config_file(config, defaults):
	up_date={}
	for category in config.keys():
		feature_list = config[category]
		for i in range(0, len(feature_list)):
			if category in defaults["general"] and feature_list[i] == "None": 
				config[category][i]=defaults["general"][category];
				print "[!] index %s: None value in general categoty %s -> updated to default value %s" %(config['index'][i], category, defaults["general"][category] )
			elif category in defaults:
				if feature_list[i] in defaults[category]:
					for element in defaults[category][feature_list[i]]:
						print "[!] index %s: updating config file with specific prmts %s for category %s::%s " %(config['index'][i], element, category, feature_list[i])
						if element in up_date: up_date[element].append(defaults[category][feature_list[i]][element])
						else: up_date[element]=[]; up_date[element].append(defaults[category][feature_list[i]][element])
				elif feature_list[i] == "None": config[category]=defaults[category]
	config.update(up_date)
	print "[$] default config file loaded"
	return(config)

user_config_dict=load_user_config_file(user_json_configuration_file)
default_configuration_dict=load_default_config_file(default_configuration_file)
update_config_dict=update_user_config_file(user_config_dict, default_configuration_dict)

#print final_config_file dictionary for debugging
tsv_output="#debugging\n#category\tfeaturesList\n"; final_config_file=update_config_dict
for element in final_config_file: tsv_output+="%s\t%s\n" %(element, final_config_file[element])
print tsv_output, "**done**"

## NEW PART FOR EXPORT OF RESULTS
json_out = json.dumps(update_config_dict)
ff = open("updated_user_configuration_file.json","w+")
ff.write(json_out)
ff.close()
