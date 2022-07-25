Files associated with python notebooks:

01_unmatched_foodcodes_fndds_asa_export.ipynb

  Input:
	
    fl100_recalls_qcd.csv - ASA24 diet recalls passing quality control
    fndds2018.csv - Food and nutrition database for dietary studies 2017-2018
		
  Output:
	
    fndds_missing_foods.csv - Foods that were missing from FNDDS 2017-2018 after mapping from versions 4.1 and 2011-2012 which were used with ASA24 v 2014 and 2016 during original data collection.

02_ingredientize_mixed_foods.ipynb

  Input:
	
    fl100_recalls_qcd.csv - ASA24 diet recalls passing quality control
    fndds_missing_foods_replaced.csv - Manually updated foodcodes linking discontinued foodcodes from v 4.1 and 2011-2012 to v 2017-2018
		
  Output:
	
    asa_fndds_matched_nndc_120721.csv - Ingredientized dataset from merging asa24 recalls with fndds2018
  
03_ingredient_code_remap

  Input:
	
    fndds2018.csv - Food and nutrition database for dietary studies 2017-2018
    asa_fndds_matched_nndc_120721.csv - Ingredientized dataset from merging asa24 recalls with fndds2018
  Output:
	
    asa_recode.csv through asa_recode6.csv - interative remapping of ingredient codes from 8-digit foodcodes
    ingred_code_remapped_102021.csv - Manually combined dataset from 03_ingredient_code_remap.ipynb outputs asa_recode.csv through asa_recode6.csv
    
04_map_ingredient_nutrient_values.ipynb

  Input:
	
    asa_fndds_matched_nndc_120721.csv -  Ingredientized dataset from merging asa24 recalls with fndds2018
    ingred_code_remapped_102021.csv - Manually combined dataset from 03_ingredient_code_remap.ipynb outputs asa_recode.csv through asa_recode6.csv
    fndds_2018_ingredient_carbohydrate_values.csv - FNDDS ingredients nutrient database. Provides carbohydrate fiber and energy (kcal) values to assess nutrient intakes from ingredientized data
		
  Output:
	
    ingredient_fiber_carb_weights_nndc_120721.csv - Ingredientized asa24 food recalls with nutrient values for carbs, fiber and energy (kcal)
    
06_merge_asa24_glycopedia_all_matches.ipynb

  Input:
	
    ingredient_fiber_carb_weights_nndc_120721.csv - Ingredientized asa24 food recalls with nutrient values for carbs, fiber and energy (kcal)
    asa_glycan_foods_subject_filtered_qcd_122221.csv - Manually matched foods/ingredients from asa24 food recalls to glycopedia
    glycopedia_wet_wt_040722.csv - Monosaccharide values (wet wt. g/g) for foods/ingredients
		
  Output:
	
    all_items_cal_adjusted_041122.csv - Intakes of monosaccharides energy adjusted (per 1000 kcal)
    all_items_unadjusted_041122.csv - Intakes of monosaccharides (not energy adjusted)
    
