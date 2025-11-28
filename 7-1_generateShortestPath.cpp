#include <Rcpp.h>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
using namespace Rcpp;


// [[Rcpp::export]]
CharacterVector shortest_path(List x) {
	//Construct unordered set for output
	std::unordered_set<std::string> output_strset;
	//Construct vector of empty strings to fill them with variants
	std::vector<std::string> variant_strvec;
	//Construct strings to store the start and end
	std::string start_str = "";
	std::string end_str = "";
	//Loop through the input list to get coordinates of the starting point
	for (int i = 0; i < x.length(); ++i) {
			IntegerVector start_ivec = x[i];
		//Loop through the input list to get coordinates of the ending point
		for (int j = 0; j < x.length(); ++j) {
			IntegerVector end_ivec = x[j];
			//Loop through the integerVectors, compare bits and memoize unequal positions
			IntegerVector positions;
			for (int k = 0; k < start_ivec.length(); ++k) {
				if (start_ivec[k] != end_ivec[k]) {
					positions.push_back(k);
					variant_strvec.push_back( "" );
				}
			}
			//Loop through the non-matching positions
			for (int l = 0; l < positions.length(); ++l) {
				//Loop through coordinates of the starting point and save all the variants
				for (int m = 0; m < start_ivec.length(); ++m) {
					//Consider this variant
					if (positions[l] == m) {
						start_ivec[m] = end_ivec[m];
					}
					variant_strvec[l] = variant_strvec[l] + std::to_string(start_ivec[m]);
					//Add dlm
					if (m < start_ivec.length() - 1) {
						variant_strvec[l] = variant_strvec[l] + "-";
					}
				}
			}
			//Make current results unique
			output_strset.insert(variant_strvec.begin(), variant_strvec.end());
			//Clear
			variant_strvec.clear();
			end_str = "";
			start_str = "";
		}
	}
	return Rcpp::wrap(output_strset);
}