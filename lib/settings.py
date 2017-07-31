from decimal import Decimal
#--------------------------------------------------------------------------------------------------#
#                               DEFAULT WEIGHTS                                                    #
#--------------------------------------------------------------------------------------------------#
weights = {
	'gap':'0.99',
	'overlap':'0.95',
	'switch':'0.05',
	'min_orf_length':'90'
	}

start_weight = {'ATG':Decimal('1'), 'CAT':Decimal('1'),
		'GTG':Decimal('0.5'), 'CAC':Decimal('0.5'), 
		'TTG':Decimal('0.05'), 'CAA':Decimal('0.05')}
