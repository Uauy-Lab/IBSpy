import sys
import IBSpy
from IBSpy.IBSpy_options import get_options


def main():
	options = get_options()
	prefs = options.preferences
	for k, v in prefs.items():
		print(f"{k}:\t{v}", file=sys.stderr)
	ibspy_results = IBSpy.IBSpyResultsSet(options= options)
	mv = ibspy_results.values_matrix.merged_values
	#ibspy_results.values_matrix.save()
if __name__ == '__main__':
	main()