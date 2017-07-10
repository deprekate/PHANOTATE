from decimal import Decimal
#--------------------------------------------------------------------------------------------------#
#                               DEFAULT WEIGHTS                                                    #
#--------------------------------------------------------------------------------------------------#
weights = {
	'gap':'0.99',
	'overlap':'0.95',
	'switch':'0.01',
	'min_orf_length':'70'
	}

start_weight = {'ATG':Decimal('5'), 'CAT':Decimal('5'),
		'GTG':Decimal('3'), 'CAC':Decimal('3'), 
		'TTG':Decimal('0.5'), 'CAA':Decimal('0.5')}
