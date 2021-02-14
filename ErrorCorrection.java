package Aligner;

import java.util.ArrayList;

public class ErrorCorrection {
	
	public int sub=-1, gap=-1, match=1, globalmax=0;
	char[] alphabet;
	byte[][] substituted_strings, chars_per_pointer;
	byte[] current_chars;
	int beta, delta;
	int[] factorials, lengths, scores, pointers;
	double subchance, gapchance;
	String str1, str2, result;
	ArrayList<double[]> dict_values;
	ArrayList<ArrayList<String>> inserts;
	
	ErrorCorrection (ArrayList<String> reads,int index, char[] alphabet_,int cov_estimate,double subchance_,double gapchance_) {
		subchance = subchance_;
		gapchance = gapchance_;
		str1 = ","+reads.get(index);
		dict_values = new ArrayList<double[]>();
		inserts = new ArrayList<ArrayList<String>>();
		for (int i=0; i<reads.get(index).length(); i++) {
			dict_values.add(new double[alphabet_.length+1]); // index : 0 is DEL, last is INS
			inserts.add(new ArrayList<String>());
		}
		result = "";
		alphabet = alphabet_;
		// score read[index] against every other read .. save mapped symbol/DEL/INS at each position, weighed with score
		for (int n=0; n<reads.size(); n++) {
			globalmax = 0;
			if (n==index) continue;
			str2 = ","+reads.get(n);
			prepare();
			for (int i=0; i<delta; i++) score(i);
			int j=0;
			int gaps=0, subs=0, length=0;
			for (j=globalmax; j!=pointers[j]; j=pointers[j],length++) {
				if (j-pointers[j] == 1 || j-pointers[j] == lengths[0]) gaps++;
				if (j-pointers[j] == lengths[0]+1 && scores[j] == scores[pointers[j]]-1) subs++;
			}
			double d = Math.exp(weight(length,gaps,subs)) * scores[globalmax]  		// <- is this an error? (e^weight) * score
					* scores[globalmax]/Math.min(str1.length()-1,str2.length()-1); 	// 					 or e^(weight * score)
			
			//*(length-(str2.length()-length)));
			for (int i=globalmax; i!=pointers[i]; i=pointers[i]) {
				int r = i%lengths[0], c = i/lengths[0];
				if (i-pointers[i] == lengths[0]+1) {
					dict_values.get(r-1)[substituted_strings[1][c]] += d;
					// System.out.println("Sub or match");
				}
				else if (i-pointers[i] == 1) {
					dict_values.get(r-1)[0] += d;
					// System.out.println("DEL");
				}
				else if (i-pointers[i] == lengths[0]) { 	// INS in str1
					// System.out.println("INS");
					String INS = "";
					for (;(i-pointers[i]) == lengths[0]; i=pointers[i]) 
						INS = str2.charAt(i/lengths[0]) + INS;
					inserts.get(r-1).add(INS);
					dict_values.get(r-1)[alphabet_.length] += d;
				}
			}
		}
		
		// find most likely inserts from saved differences
		for (int i=inserts.size()-1; i>=0; i--) {
			double[] tmp = new double[alphabet_.length+1]; 	// copy i-th entry of dict_values-list minus INS score
			int symb = -1; 		// index of most likely symbol/DEL/INS
			double maxc = 0; 	// tmp for assessing most likely symbol
			
			for (int k=0,limit=tmp.length-1; k<=limit; k++) {
				if (k!=limit) tmp[k] = dict_values.get(i)[k]-dict_values.get(i)[limit];
				else {
					tmp[k] = dict_values.get(i)[k];
					dict_values.get(i)[k] = 0;
				}
				if (tmp[k] > maxc) {
					maxc = tmp[k];
					symb = k;
				}
			}
			
			if (symb == tmp.length-1) { 		// if INS
				dict_values.add(i+1,tmp.clone());
				inserts.add(i+1,new ArrayList<String>());
			}
			
			String INS = "";
			for (int k=0; ;k++) {
				int[] symbols = new int[alphabet_.length];
				int maxcounter = 0;
				char maxsymb = '\0';
				for (int j=0; j<inserts.get(i).size(); j++) {
					char insrt = '\0';
					int m = 0;
					if (k<inserts.get(i).get(j).length()) {
						m = index_in_alphabet(inserts.get(i).get(j).charAt(k));
						symbols[m]++;
						insrt = inserts.get(i).get(j).charAt(k);
					} else {
						symbols[0]++;
						insrt = '\0';
					}
					if (maxcounter<=symbols[m]) {
						maxcounter = symbols[m];
						maxsymb = insrt;
					}
				}
				if (maxsymb!='\0') INS += maxsymb;
				else break;
				
			}
			if (symb == tmp.length-1) { 			// if INS
				inserts.get(i+1).add(0,INS);
			}
			
		}
		
		for (int k=0; k<dict_values.size(); k++) {
			int symb = -1;
			double maxc = 0;
			for (int j=0; j<=alphabet_.length; j++) {
				if (dict_values.get(k)[j]>maxc) {
					maxc = dict_values.get(k)[j];
					symb = j;
				}
			}
			if (symb == alphabet_.length) result += inserts.get(k).get(0);
			else if (symb > 0) result += alphabet[symb];
		}
		
		reads.set(index,result);
	}
	
	// only applicable for initial alphabet (without the comma)
	public int index_in_alphabet(char c) {
		for (int i=0; i<alphabet.length; i++)
			if (alphabet[i]==c) return i-1;
		return -1;
	}
	
	private double weight(int length, int gaps, int subs) {
		double gapweight = 1, subweight = 1;
		for (int i=0; i<=gaps; i++) gapweight -= 1. * binomial_coeff(length,i)*Math.pow(gapchance,i)*Math.pow(1-gapchance,length-i);
		for (int i=0; i<=subs; i++) gapweight -= 1. * binomial_coeff(length,i)*Math.pow(subchance,i)*Math.pow(1-subchance,length-i);
		return subweight*gapweight;
	}
	
	private long binomial_coeff(int n, int k) {
		long res = 1;
		for (int i=n-k+1; i<=n; i++) res *= i;
		for (int i=1; i<=k; i++) res /= i;
		return res;
	}
	
	public void prepare() {
		
		// generate number arrays from input strings 
		substituted_strings = new byte[2][];
		substituted_strings[0] = new byte[str1.length()];
		for (int j=0; j<str1.length(); j++)
			for (byte n=0; n<alphabet.length; n++)
				if (str1.charAt(j) == alphabet[n])
					substituted_strings[0][j] = n;
		substituted_strings[1] = new byte[str2.length()];
		for (int j=0; j<str2.length(); j++)
			for (byte n=0; n<alphabet.length; n++)
				if (str2.charAt(j) == alphabet[n])
					substituted_strings[1][j] = n;
		
		beta = alphabet.length;
		chars_per_pointer = new byte[3][beta];
		factorials = new int[2];
		factorials[0] = 1;
		factorials[1] = substituted_strings[0].length;
		lengths = new int[2];
		for (int k=0; k<2; lengths[k]=substituted_strings[k].length,k++);
		delta = factorials[1]*lengths[1];
		scores = new int[delta];
		pointers = new int[delta];
		current_chars = new byte[2];
		
	}
	
	public void score(int index) {
		chars_per_pointer[0] = new byte[beta];
		int epsilon = 0;	// maxpointer
		
		for (int i=0,k=index,j=0; i<2; k/=lengths[i],i++) {
			j = k%lengths[i];
			if (j>0) {
				epsilon += Math.pow(2,i);
				current_chars[i] = substituted_strings[i][j];
				chars_per_pointer[0][substituted_strings[i][j]]++;
			} else current_chars[i] = -1;
		}
		
		int[] current_pointers = new int[3];
		int[] pointer_scores = new int[3];
		current_pointers[0] = epsilon;
		
		// get pointer score for all pointers except maxpointer
		double zeta = 0; 	// zeros in maxpointer
		for (int i=1,k=1,j=0; i<4; i<<=1,j++) {
			if ((i&epsilon)>0) {
				int n = k;
				while (n>0) {
					n--;
					int m = current_pointers[n] - i;
					if (m>0) {
						current_pointers[k] = m;
						int q = 0;
						for (int p=0; p<beta; p++) {
							byte t = chars_per_pointer[n][p];
							if (p==current_chars[j]) t--;
							q += (t*(t-1))/2;
							chars_per_pointer[k][p] = t;
						}
						pointer_scores[k] = pointer_scores[n] + gap + (match-sub) * q;
						k++;
					}
				}
			} else zeta += gap;
		}
		
		// get score for maxpointer
		int eta = 0;
		for (int i=0,k; i<beta; k=chars_per_pointer[0][i],eta+=(k*(k-1))/2,i++);
		pointer_scores[0] = (match-sub) * eta;
		
		for (int i=0; i<3; i++) {
			int j = 0;
			for (int k=1,n=0; k<4; n++,k<<=1)
				if ((k&current_pointers[i])>0) j += factorials[n];
			current_pointers[i] = index - j;
			int charsum = 0;
			for (int k=0; k<chars_per_pointer[i].length; charsum+=chars_per_pointer[i][k],k++);
			pointer_scores[i] += zeta + sub * (charsum*(charsum-1))/2;
		}
		
		if (index>0) {
			int max = scores[current_pointers[0]] + pointer_scores[0];
			int max_i = 0;
			for (int i=1; i<3 && current_pointers[i]!=index; i++) 
				if (scores[current_pointers[i]] + pointer_scores[i] > max) {
					max = scores[current_pointers[i]] + pointer_scores[i];
					max_i = i;
				}
			if (index%lengths[0]>0 && index/lengths[0]>0) {
				scores[index] = max;
				pointers[index] = current_pointers[max_i];
				if (max>=scores[globalmax] && (index%lengths[0]+1 == str1.length() || index/lengths[0]+1 == str2.length()))
					globalmax = index;
			} else {
				scores[index] = 0;
				pointers[index] = index;
			}
		} else {
			scores[0] = 0;
			pointers[0] = 0;
		}
	}
	
	/*
	ErrorCorrection (ArrayList<String> reads,int index, char[] alphabet_,int cov_estimate,double subchance_,double gapchance_) {
		subchance = subchance_;
		gapchance = gapchance_;
		str1 = ","+reads.get(index);
		dict_values = new ArrayList<double[]>();
		inserts = new ArrayList<ArrayList<String>>();
		for (int i=0; i<reads.get(index).length(); i++) {
			dict_values.add(new double[alphabet_.length+1]); // index : 0 is DEL, last is INS
			inserts.add(new ArrayList<String>());
		}
		result = "";
		alphabet = alphabet_;
		// score read[index] against every other read .. save mapped symbol/DEL/INS at each position, weighed with score
		for (int n=0; n<reads.size(); n++) {
			globalmax = 0;
			if (n==index) continue;
			str2 = ","+reads.get(n);
			prepare();
			for (int i=0; i<delta; i++) score(i);
			int j=0;
			int gaps=0, subs=0, length=0, i=0;
			for (j=globalmax; j!=pointers[j]; j=pointers[j],length++) {
				if (j-pointers[j] == 1 || j-pointers[j] == lengths[0]) gaps++;
				if (j-pointers[j] == lengths[0]+1 && scores[j] == scores[pointers[j]]-1) subs++;
			}
			double d = Math.exp(scores[globalmax]
								- Math.min(j/lengths[0],j%lengths[0])
								- Math.min(str2.length()-1-globalmax/lengths[0],
										   str1.length()-1-globalmax%lengths[0]));
			
			//*(length-(str2.length()-length)));
			for (int i=globalmax; i!=pointers[i]; i=pointers[i]) {
				int r = i%lengths[0], c = i/lengths[0];
				if (i-pointers[i] == lengths[0]+1) {
					dict_values.get(r-1)[substituted_strings[1][c]] += d;
					// System.out.println("Sub or match");
				}
				else if (i-pointers[i] == 1) {
					dict_values.get(r-1)[0] += d;
					// System.out.println("DEL");
				}
				else if (i-pointers[i] == lengths[0]) { 	// INS in str1
					// System.out.println("INS");
					String INS = "";
					for (;(i-pointers[i]) == lengths[0]; i=pointers[i]) 
						INS = str2.charAt(i/lengths[0]) + INS;
					inserts.get(r-1).add(INS);
					dict_values.get(r-1)[alphabet_.length] += d;
				}
			}
		}
		
		// find most likely inserts from saved differences
		for (int i=inserts.size()-1; i>=0; i--) {
			double[] tmp = new double[alphabet_.length+1]; 	// copy i-th entry of dict_values-list minus INS score
			int symb = -1; 		// index of most likely symbol/DEL/INS
			double maxc = 0; 	// tmp for assessing most likely symbol
			
			for (int k=0,limit=tmp.length-1; k<=limit; k++) {
				if (k!=limit) tmp[k] = dict_values.get(i)[k]-dict_values.get(i)[limit];
				else {
					tmp[k] = dict_values.get(i)[k];
					dict_values.get(i)[k] = 0;
				}
				if (tmp[k] > maxc) {
					maxc = tmp[k];
					symb = k;
				}
			}
			
			if (symb == tmp.length-1) { 		// if INS
				dict_values.add(i+1,tmp.clone());
				inserts.add(i+1,new ArrayList<String>());
			}
			
			String INS = "";
			for (int k=0; ;k++) {
				int[] symbols = new int[alphabet_.length];
				int maxcounter = 0;
				char maxsymb = '\0';
				for (int j=0; j<inserts.get(i).size(); j++) {
					char insrt = '\0';
					int m = 0;
					if (k<inserts.get(i).get(j).length()) {
						m = index_in_alphabet(inserts.get(i).get(j).charAt(k));
						symbols[m]++;
						insrt = inserts.get(i).get(j).charAt(k);
					} else {
						symbols[0]++;
						insrt = '\0';
					}
					if (maxcounter<=symbols[m]) {
						maxcounter = symbols[m];
						maxsymb = insrt;
					}
				}
				if (maxsymb!='\0') INS += maxsymb;
				else break;
				
			}
			if (symb == tmp.length-1) { 			// if INS
				inserts.get(i+1).add(0,INS);
			}
			
		}
		
		for (int k=0; k<dict_values.size(); k++) {
			int symb = -1;
			double maxc = 0;
			for (int j=0; j<=alphabet_.length; j++) {
				if (dict_values.get(k)[j]>maxc) {
					maxc = dict_values.get(k)[j];
					symb = j;
				}
			}
			if (symb == alphabet_.length) result += inserts.get(k).get(0);
			else if (symb > 0) result += alphabet[symb];
		}
		
		reads.set(index,result);
	}
	
	// only applicable for initial alphabet (without the comma)
	public int index_in_alphabet(char c) {
		for (int i=0; i<alphabet.length; i++)
			if (alphabet[i]==c) return i-1;
		return -1;
	}
	
	private double weight(int length, int gaps, int subs) {
		double gapweight = 1. * binomial_coeff(length,gaps)*Math.pow(subchance+gapchance,gaps)*Math.pow(1-subchance+gapchance,length-gaps);
		double subweight = 1. * binomial_coeff(length,subs)*Math.pow(subchance+gapchance,subs)*Math.pow(1-subchance+gapchance,length-subs);
		return subweight*gapweight;
	}
	
	private long binomial_coeff(int n, int k) {
		long res = 1;
		for (int i=n-k+1; i<=n; i++) res *= i;
		for (int i=1; i<=k; i++) res /= i;
		return res;
	}
	
	public void prepare() {
		
		// generate number arrays from input strings 
		substituted_strings = new byte[2][];
		substituted_strings[0] = new byte[str1.length()];
		for (int j=0; j<str1.length(); j++)
			for (byte n=0; n<alphabet.length; n++)
				if (str1.charAt(j) == alphabet[n])
					substituted_strings[0][j] = n;
		substituted_strings[1] = new byte[str2.length()];
		for (int j=0; j<str2.length(); j++)
			for (byte n=0; n<alphabet.length; n++)
				if (str2.charAt(j) == alphabet[n])
					substituted_strings[1][j] = n;
		
		beta = alphabet.length;
		chars_per_pointer = new byte[3][beta];
		factorials = new int[2];
		factorials[0] = 1;
		factorials[1] = substituted_strings[0].length;
		lengths = new int[2];
		for (int k=0; k<2; lengths[k]=substituted_strings[k].length,k++);
		delta = factorials[1]*lengths[1];
		scores = new int[delta];
		pointers = new int[delta];
		current_chars = new byte[2];
		
	}
	
	public void score(int index) {
		chars_per_pointer[0] = new byte[beta];
		int epsilon = 0;	// maxpointer
		
		for (int i=0,k=index,j=0; i<2; k/=lengths[i],i++) {
			j = k%lengths[i];
			if (j>0) {
				epsilon += Math.pow(2,i);
				current_chars[i] = substituted_strings[i][j];
				chars_per_pointer[0][substituted_strings[i][j]]++;
			} else current_chars[i] = -1;
			if (substituted_strings[i][j]==0) {
				scores[index] = 0;
				pointers[index] = index;
				return;
			}
		}
		
		int[] current_pointers = new int[3];
		int[] pointer_scores = new int[3];
		current_pointers[0] = epsilon;
		
		// get pointer score for all pointers except maxpointer
		double zeta = 0; 	// zeros in maxpointer
		for (int i=1,k=1,j=0; i<4; i<<=1,j++) {
			if ((i&epsilon)>0) {
				int n = k;
				while (n>0) {
					n--;
					int m = current_pointers[n] - i;
					if (m>0) {
						current_pointers[k] = m;
						int q = 0;
						for (int p=0; p<beta; p++) {
							byte t = chars_per_pointer[n][p];
							if (p==current_chars[j]) t--;
							q += (t*(t-1))/2;
							chars_per_pointer[k][p] = t;
						}
						pointer_scores[k] = pointer_scores[n] + gap + (match-sub) * q;
						k++;
					}
				}
			} else zeta += gap;
		}
		
		// get score for maxpointer
		int eta = 0;
		for (int i=0,k; i<beta; k=chars_per_pointer[0][i],eta+=(k*(k-1))/2,i++);
		pointer_scores[0] = (match-sub) * eta;
		
		for (int i=0; i<3; i++) {
			int j = 0;
			for (int k=1,n=0; k<4; n++,k<<=1)
				if ((k&current_pointers[i])>0) j += factorials[n];
			current_pointers[i] = index - j;
			int charsum = 0;
			for (int k=0; k<chars_per_pointer[i].length; charsum+=chars_per_pointer[i][k],k++);
			pointer_scores[i] += zeta + sub * (charsum*(charsum-1))/2;
		}
		
		if (index>0) {
			int max = scores[current_pointers[0]] + pointer_scores[0];
			int max_i = 0;
			for (int i=1; i<3 && current_pointers[i]!=index; i++) 
				if (scores[current_pointers[i]] + pointer_scores[i] > max) {
					max = scores[current_pointers[i]] + pointer_scores[i];
					max_i = i;
				}
			if (max>0) {
				scores[index] = max;
				pointers[index] = current_pointers[max_i];
				if (max>=scores[globalmax]) globalmax = index;
			} else {
				scores[index] = 0;
				pointers[index] = index;
			}
		} else {
			scores[0] = 0;
			pointers[0] = 0;
		}
	}
	**/

}
