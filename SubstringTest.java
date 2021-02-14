package Aligner;

public class SubstringTest {
	
	private int sub=-1, gap=-1, match=1, globalmax=0;
	private char[] alphabet;
	private byte[][] substituted_strings, chars_per_pointer;
	private double subchance, gapchance;
	private byte[] current_chars;
	private int beta, delta;
	private int[] factorials, lengths, scores, pointers;
	private String str1, str2;
	public boolean result = false,done = false;

	SubstringTest (String st1, String st2, char[] alphabet_, double subchance_, double gapchance_) {
		// len(st1) <= len(st2) !!
		if (st2.indexOf(st1)>=0) result = true;
		if (subchance_<0.000000000000000001 && gapchance_<0.000000000000000001) {
			done = true;
			return;
		}
		else {
			str1 = ","+st1;
			str2 = ","+st2;
			alphabet = alphabet_;
			subchance = subchance_;
			gapchance = gapchance_;
			prepare();
			for (int i=0; i<delta; i++) score(i);
			result = check_if_substring();
			
		}
		done = true;
	}
	
	private boolean check_if_substring() {
		int gaps=0, subs=0, length=0, i=0;
		for (i=globalmax; i!=pointers[i]; i=pointers[i],length++) {
			if (i-pointers[i] == 1 || i-pointers[i] == lengths[0]) gaps++;
			if (i-pointers[i] == lengths[0]+1 && scores[i] == scores[pointers[i]]-1) subs++;
		}
		gaps += i%lengths[0] + str1.length()-1 - globalmax%lengths[0];
		
		if (gaps+subs>error_cutoff(length)) return false;
		return true;
	}
	
	private void prepare() {
		
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
	
	private int error_cutoff(int length) {
		double sum=0;
		int i=0;
		while (true) {
			double add = 1. * binomial_coeff(length,i)*Math.pow(subchance+gapchance,i)*Math.pow(1-subchance+gapchance,length-i);
			sum += add;
			if (add/(sum+add) < 0.001) break;
			i++;
		}
		return i;
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
	
	private void score(int index) {
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

}
