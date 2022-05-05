import sys
from webbrowser import Chrome
import IBSpy
from IBSpy.IBSpy_options import get_options

def main():
	options = get_options()
	prefs = options.preferences
	for k, v in prefs.items():
		print(f"{k}:\t{v}", file=sys.stderr)
	ibspy_results = IBSpy.IBSpyResultsSet(options= options)
	chromosome, start, end, assembly =  options.region
	
	window = ibspy_results.mapped_window_tabix(chromosome, start, end, assembly=assembly)
	# print(f"[_function_window_wrapper_tabix]About to run region: {chromosome}:{start}-{end}")
	window = window.drop_duplicates(keep="last")
	window.to_csv(options.csv_output, sep="\t",index=False)

	#ibspy_results.values_matrix.save()
if __name__ == '__main__':
	main()