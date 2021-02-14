package Aligner;

class ErrorCount {
	
	static public int sub=1,gap=1;
	static char[] alphabet;
	static byte[][] substituted_strings,chars_per_pointer;
	static byte[] current_chars;
	static int alpha,beta,gamma,delta;
	static int[] factorials,lengths;
	static double[] scores;
	static int[] pointers;
	
	static int[][] trace;
	static int[][] aligned_pointers;

	public ErrorCount(String seq1, String seq2,int line_width, boolean print_alignment) {
	
		String str1 = ","+seq1+",";
		String str2 = ","+seq2+",";

		prepare(str1,str2);
		for (int i=0; i<delta; i++) score(i);
		traceback();
		double error_rate = print(line_width,print_alignment);
		System.out.printf("\n                 Error rate : \t%.3f%s",error_rate*100,"%");

	}
	
	static public void prepare(String str1, String str2) {
		
		String[] input = {str1,str2};
		String[] input2;
		
		// determine string alphabet
		alphabet = new char[0];
		for (int k=0; k<input.length; k++)
			for (int j=0; j<input[k].length(); j++) {
				boolean in_alphabet = false;
				for (int l=0; l<alphabet.length; l++)
					if (input[k].charAt(j) == alphabet[l])
						in_alphabet = true;
				if (!in_alphabet) {
					input2 = new String[alphabet.length];
					for (int n=0; n<alphabet.length; n++)
						input2[n] = String.valueOf(alphabet[n]);
					alphabet = new char[alphabet.length+1];
					for (int n=0; n<input2.length; n++)
						alphabet[n] = input2[n].charAt(0);
					alphabet[alphabet.length-1] = input[k].charAt(j);
				}
			}
		
		// generate number arrays from input strings 
		substituted_strings = new byte[input.length][];
		for (int k=0; k<input.length; k++) {
			substituted_strings[k] = new byte[input[k].length()];
			for (int j=0; j<input[k].length(); j++)
				for (byte n=0; n<alphabet.length; n++)
					if (input[k].charAt(j) == alphabet[n])
						substituted_strings[k][j] = n;
		}
		
		// nr of strings (v13)
		alpha = input.length;
		// size of alphabet (v12)
		beta = alphabet.length;
		// max number of pointers
		gamma = (int)(Math.pow(2, alpha))-1;
		
		// l5 (char count for each pointer)
		chars_per_pointer = new byte[gamma][beta];
		// l6 ("factorials")
		factorials = new int[alpha];
		factorials[0] = 1;
		for (int k=1; k<alpha; factorials[k]=factorials[k-1]*substituted_strings[k-1].length,k++);
		
		// length of the strings
		lengths = new int[alpha];
		for (int k=0; k<alpha; lengths[k]=substituted_strings[k].length,k++);
		
		// size of scoring and backpointer matrix
		delta = factorials[alpha-1]*lengths[alpha-1];
		scores = new double[delta];
		pointers = new int[delta];
		
		// chars at current position
		current_chars = new byte[alpha];
		
	}
	
	public static void score(int index) {
		chars_per_pointer[0] = new byte[beta];
		int epsilon = 0;	// maxpointer
		
		for (int i=0,k=index,j=0; i<alpha; k/=lengths[i],i++) {
			j = k%lengths[i];
			if (j>0) {
				epsilon += Math.pow(2,i);
				current_chars[i] = substituted_strings[i][j];
				chars_per_pointer[0][substituted_strings[i][j]]++;
			} else current_chars[i] = -1;
		}
		
		int[] current_pointers = new int[gamma];
		double[] pointer_scores = new double[gamma];
		current_pointers[0] = epsilon;
		
		// get pointer score for all pointers except maxpointer
		double zeta = 0; 	// zeta = nr of zeros in maxpointer * gap
		for (int i=1,k=1,j=0; i<1<<alpha; i<<=1,j++) { 			// this loop iterates over all the strings; j corresponds to i as in : i = 2^j 
			if ((i&epsilon)>0) { 		// if the maxpointer is 1 in dimension i, create a new pointer which is 0 in this dimension from every pointer already created  
				int n = k;  			// n iterates over all pointers that were created already
				while (n>0) {
					n--;
					int m = current_pointers[n] - i; 			// subtract 2^(j-1) --> creates 0 at position j in binary number
					if (m>0) {
						current_pointers[k] = m; 				// store difference as new pointer
						int q = 0; 								// q is to be number of matching pairs
						for (int p=0; p<beta; p++) { 				// iterate p over length of alphabet
							byte t = chars_per_pointer[n][p]; 		// t := number of chars of type alphabet[p] that are traversed by pointer[n]
							if (p==current_chars[j]) t--; 			// subtract 1 from number of chars of type alphabet[j]
							q += (t*(t-1))/2; 						// q is increased by number of matching pairs for each char of the alphabet 
							chars_per_pointer[k][p] = t;
						}
						// subscore is subtracted for each match, bc later its added for every possible pair
						pointer_scores[k] = pointer_scores[n] + gap - sub * q;
						k++;
					}
				}
			} else zeta += gap;
		}
		
		// get score for maxpointer
		int eta = 0;
		for (int i=0,k; i<beta; k=chars_per_pointer[0][i],eta+=(k*(k-1))/2,i++);
		pointer_scores[0] = -sub * eta;
		
		for (int i=0; i<gamma; i++) {
			int j = 0;
			for (int k=1,n=0; k<1<<alpha; n++,k<<=1)
				if ((k&current_pointers[i])>0) j += factorials[n];
			current_pointers[i] = index - j;
			int charsum = 0;
			for (int k=0; k<chars_per_pointer[i].length; charsum+=chars_per_pointer[i][k],k++);
			pointer_scores[i] += zeta + sub * (charsum*(charsum-1))/2;
		}
		
		if (index>0) {
			double max = scores[current_pointers[0]] + pointer_scores[0];
			int max_i = 0;
			for (int i=1; i<gamma && current_pointers[i]!=index; i++) 
				if (scores[current_pointers[i]] + pointer_scores[i] < max) {
					max = scores[current_pointers[i]] + pointer_scores[i];
					max_i = i;
				}
			scores[index] = max;
			pointers[index] = current_pointers[max_i];
		} else {
			scores[0] = 0;
			pointers[0] = 0;
		}
	}
	
	public static void traceback() {
		trace = new int[1][alpha];
		int[][] tmp_trace;
		int i = delta-1,k=0;
		boolean last = false;
		
		// trace back
		do {
			if (pointers[i]==i) last = true;
			for (int j=0; j<alpha; j++)
				trace[k][j] = (i/factorials[j])%lengths[j];
			k++;
			i = pointers[i];
			tmp_trace = trace;
			trace = new int[trace.length+1][alpha];
			for (int j=0; j<tmp_trace.length; trace[j]=tmp_trace[j], j++);
		} while (!last);
		
		trace = new int[tmp_trace.length][alpha];
		// revert trace
		for (int j=tmp_trace.length-1,n=0; j>-1; j--,n++) 
			trace[j] = tmp_trace[n];
		
		// create pointers to trace
		aligned_pointers = new int[trace.length-1][alpha];
		for (int j=0; j<trace.length-1; j++) {
			int[] pointer = new int[alpha];
			for (int n=0; n<alpha; pointer[n]=trace[j+1][n]-trace[j][n],n++);
			aligned_pointers[j] = pointer;
		}
		
		tmp_trace = trace;
		trace = new int[tmp_trace.length-2][alpha];
		for (int j=0; j<trace.length; trace[j]=tmp_trace[j+1],j++);
		tmp_trace = aligned_pointers;
		aligned_pointers = new int[trace.length][alpha];
		for (int j=0; j<tmp_trace.length-1; aligned_pointers[j]=tmp_trace[j],j++);
		
	}
	
	// returns error rate
	public static double print(int line_width, boolean print_alignment) {
		String[][] aligned_strings = new String[(alpha*(alpha-1))/2][3];
		for (int i=0,n=0; i<alpha-1; i++)
			for (int k=i+1; k<alpha; k++,n++) {
				aligned_strings[n][0] = "";
				aligned_strings[n][1] = "";
				aligned_strings[n][2] = "";
				for (int j=0; j<trace.length; j++) {
					char a,b,c = 0;
					if (aligned_pointers[j][i] == 1) a = alphabet[substituted_strings[i][trace[j][i]]];
					else a = '_';
					if (aligned_pointers[j][k] == 1) b = alphabet[substituted_strings[k][trace[j][k]]];
					else b = '_';
					boolean append = true;
					if (a==b) {
						if (a=='_') append = false;
						else c = '|';
					} else c = '*';
					if (append) {
						aligned_strings[n][0] += String.valueOf(a);
						aligned_strings[n][1] += String.valueOf(c);					
						aligned_strings[n][2] += String.valueOf(b);
					}
				}
			}
		
		double errors = 0;
		for (int k=0; k<aligned_strings[0][1].length(); k++) 
			if (aligned_strings[0][1].charAt(k)=='*') errors++;
		
		if (print_alignment) System.out.println("\n");
		for (String[] alstr: aligned_strings) {
			
			for (int i=0; i<(double)alstr[0].length()/line_width; i++) {
				String a=alstr[0].substring(line_width*i,(int)Math.min(line_width*(i+1),alstr[0].length()));
				String b=alstr[1].substring(line_width*i,(int)Math.min(line_width*(i+1),alstr[1].length()));
				String c=alstr[2].substring(line_width*i,(int)Math.min(line_width*(i+1),alstr[1].length()));
				String indent1="           ", indent2="           ", indent3="           ";
				if (i==0) {
					indent1 = "truth      ";
					indent2 = "           ";
					indent3 = "assembly   ";
				}
				
				if (print_alignment) System.out.println(indent1+a+System.lineSeparator()+indent2+b+System.lineSeparator()+indent3+c+System.lineSeparator()+"           ");
			}
		}
		return (errors/aligned_strings[0][0].length());
		
	}

}
