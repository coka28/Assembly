package Aligner;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Assembly {
	
	static public volatile double[][] scores,scores2,scores3,in_scores,out_scores,order;
	static public double[] rowSums,colSums;
	static public volatile String[][] merge_results;
	static String truth = "";
	static ArrayList<Character> alphabet = new ArrayList<Character>();
	static char[] alphabet_array;
	static ArrayList<String> reads = new ArrayList<String>();
	static double subchance=1e-9,
				  gapchance=1e-9,
				  sub,gap,
				  confidence_threshold,
				  speedup,
				  avrg_length = 0;
	static int threadNr=4,alphabet_size,n,m,coverage_estimate;
	static boolean print_merges = true,
				   print_init = true,
				   print_alignment = true,
				   print_progbar = true,
				   force_result = false,
				   reduce_errors=true,
				   correct_all = false;
	static public volatile double progress=0;
	static public volatile int percent=1;
	static int line_width;

	public static void main(String[] args) throws InterruptedException {
		// args to variables :
		System.setErr(System.out);
		subchance = Double.parseDouble(args[0]);
		if (subchance<=0) subchance = 10*Double.MIN_VALUE;
		if (subchance>0.1) subchance = 0.1;
		gapchance = Double.parseDouble(args[1]);
		if (gapchance<=0) gapchance = 10*Double.MIN_VALUE;
		if (gapchance>0.1) gapchance = 0.1;
		threadNr = Integer.parseInt(args[2]);
		if (args[3].equals("0")) print_merges = false;
		if (args[4].equals("0")) print_init = false;
		if (args[5].equals("1")) force_result = true;
		if (args[6].equals("0")) print_alignment = false;
		if (args[7].equals("0")) print_progbar = false;
		int error_correction_level = Integer.parseInt(args[8]);
		if (args[9].equals("1")) correct_all = true;
		confidence_threshold = Double.parseDouble(args[10]);
		line_width = Integer.parseInt(args[11]);
		sub = -Math.log(subchance);
		gap = -Math.log(gapchance);
		File read_file = new File("reads.txt");
		File truth_file = new File("truth.txt");
		try {
			BufferedReader br = new BufferedReader(new FileReader(truth_file));
			String tmp_char;
			while ((tmp_char = br.readLine()) != null) truth += tmp_char;
			br.close();
		} catch (FileNotFoundException e) {
			System.out.println("FATAL ERROR : Could not find the truth.txt file in the working directory!");
			System.exit(1);
		} catch (IOException e) {
			System.out.println("FATAL ERROR : Could not read the truth.txt file in the working directory!");
			System.exit(1);
		}
		String read_string = "";
		try {
			BufferedReader br = new BufferedReader(new FileReader(read_file));
			String tmp_char;
			read_string = "";
			while ((tmp_char = br.readLine()) != null) read_string += tmp_char;
			br.close();
		} catch (FileNotFoundException e) {
			System.out.println("FATAL ERROR : Could not open the reads.txt file in the working directory!");
			System.exit(1);
		} catch (IOException e) {
			System.out.println("FATAL ERROR : Could not read the reads.txt file in the working directory!");
			System.exit(1);
		}
		for (String read : read_string.split(",")) reads.add(read);
		get_alphabet_size();
		ArrayList<Integer> del_list = new ArrayList<Integer>();
		
		if (correct_all) {
			while (error_correction_level>0) {
				progress = 0;
				System.out.printf("Running error correction on %d reads...\n\n",reads.size());
				m = reads.size();
				ScoringQueue[] thrds = new ScoringQueue[threadNr];
				int[][][] workloads = new int[threadNr][][];
				for (int i=0,q=reads.size(); i<threadNr; i++) {
					int low=q*i/threadNr, up=q*(i+1)/threadNr;
					workloads[i] = new int[1][up-low];
					for (int k=low,j=0; k<up; j++,k++) 
						workloads[i][0][j] = k;
					thrds[i] = new ScoringQueue(workloads[i],false,false,false,i+1);
					thrds[i].start();
				}
				for (int i=0; i<threadNr; i++) thrds[i].thrd.join(0);
				error_correction_level--;
			}
		}
		
		
		for (int p=0; p<error_correction_level; p++) {
			progress = 0;
			ArrayList<String> old_reads = new ArrayList<String>(reads);
			del_list = find_substring_reads();
			coverage_estimate = (del_list.size()*2)/(reads.size()-del_list.size());
			System.out.printf("Estimated minimum coverage of %d\n",coverage_estimate);
			ArrayList<Integer> valid_indices = new ArrayList<Integer>();
			for (int i=0; i<reads.size(); i++) if (!del_list.contains(i)) valid_indices.add(i);
			System.out.printf("Running error correction on %d reads...\n\n",valid_indices.size());
			ScoringQueue[] thrds = new ScoringQueue[threadNr];
			int[][][] workloads = new int[threadNr][][];
			m = valid_indices.size();
			for (int i=0,q=valid_indices.size(); i<threadNr; i++) {
				int low=q*i/threadNr, up=q*(i+1)/threadNr;
				workloads[i] = new int[1][up-low];
				for (int k=low,j=0; k<up; j++,k++) 
					workloads[i][0][j] = valid_indices.get(k);
				thrds[i] = new ScoringQueue(workloads[i],false,false,false,i+1);
				thrds[i].start();
			}
			for (int i=0; i<threadNr; i++) thrds[i].thrd.join(0);
			if (reads.containsAll(old_reads)) break;
		}
		delete_reads(del_list);
		delete_substring_reads();
		
		sub /= Math.log(alphabet_size);
		gap /= Math.log(alphabet_size);
		n = reads.size();
		long n_ = n;
		for (String r : reads) avrg_length += r.length();
		avrg_length /= reads.size();
		double estRuntime = Math.pow(avrg_length,1.9)*2.46e-7*(n_*(n_-1)/2+n_*(n_-1)*(n_-2))/threadNr+3;
		int estHrs  = (int)(estRuntime/3600);
		int estMins = (int)(estRuntime/60)%60;
		int estSecs = (int)(estRuntime)%60;
		if (print_init) {
			System.out.printf("\n\n             Number of reads: \t%d",n);
			System.out.printf("\n         Average read length: \t%.3f\n\n",avrg_length);
		}
		System.out.print("   Estimated minimum runtime: \t");
		System.out.printf("%02d:%02d:%02d\n",estHrs,estMins,estSecs);
		try {start();}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		
	}
	
	public static void get_alphabet_size() {
		for (String r : reads)
			for (int i=0; i<r.length(); i++) 
				if (!alphabet.contains(r.charAt(i))) alphabet.add(r.charAt(i));
		alphabet_size = alphabet.size();
		alphabet_array = new char[alphabet_size+1];
		alphabet_array[0] = ',';
		for (int i=1; i<alphabet_size+1; i++) alphabet_array[i] = alphabet.get(i-1);
	}
	
	public static ArrayList<Integer> find_substring_reads() {
		ArrayList<Integer> del_list = new ArrayList<Integer>();
		for (int i=0; i<reads.size(); i++) 
			for (int k=0; k<reads.size(); k++) {
				if (i==k) continue;
				int short_index = i, long_index = k;
				if (reads.get(k).length()<reads.get(i).length()) {
					short_index = k;
					long_index = i;
				}
				if (del_list.contains(short_index)) continue;
				SubstringTest test = new SubstringTest(reads.get(short_index),reads.get(long_index),alphabet_array,subchance,gapchance);
				while (!test.done) ;
				if (test.result) {
					if (del_list.contains(long_index)) continue;
					int j=0; 
					del_list.add(j,short_index);
					break;
				}
			}
		for (int i=0; i<del_list.size()-1; i++)
			for (int k=0; k<del_list.size(); k++) 
				if (del_list.get(i)<del_list.get(i+1)) {
					int j = del_list.get(i+1);
					del_list.set(i+1,del_list.get(i));
					del_list.set(i,j);
				}
		return del_list;
	}
	
	public static void delete_substring_reads() {
		ArrayList<Integer> del = find_substring_reads();
		delete_reads(del);
	}

	public static void delete_reads(ArrayList<Integer> del_list) {
		ArrayList<String> tmp = new ArrayList<String>();
		for (int i=0; i<reads.size(); i++)
			if (!del_list.contains(i)) tmp.add(reads.get(i));
		reads = tmp;
	}
	
	public static void start() throws InterruptedException {
		long t0 = System.currentTimeMillis();
		boolean abort = false;
		n = reads.size();
		
		while (n>2) {
			progress = 0;
			percent = 1;
			System.out.println();
			scores = new double[n+1][n+1]; 		// scores from simple pairwise scoring
			scores2 = new double[n][n]; 		// rating as to which merges have higher scores than the individual reads
			scores3 = new double[n+1][n+1]; 	// rating as to which merges have other reads fitting between them and which don't
			in_scores = new double[n+1][n+1];	// relative scores in columns
			out_scores = new double[n+1][n+1];	// relative scores in rows
			merge_results = new String[n][n]; 	// matrix with the strings resulting from merging two reads
			
			ScoringQueue[] thrds = new ScoringQueue[threadNr];
			
			
			int[][][] workloads = new int[threadNr][][];
					
			for (int i=0; i<threadNr; i++) {
				int low=n*(n-1)/2*i/threadNr, up=n*(n-1)/2*(i+1)/threadNr;
				workloads[i] = new int[up-low][2];
				for (int k=low+1,j=0; k<up+1; k++,j++) {
					workloads[i][j][0] = (int) Math.ceil(-0.5+Math.sqrt(0.25+2*k));
					workloads[i][j][1] = (int) (k-workloads[i][j][0]*(workloads[i][j][0]-1)/2)-1;
				}
				thrds[i] = new ScoringQueue(workloads[i],true,false,false,i+1);
				thrds[i].start();
			}
			
			for (int i=0; i<threadNr; i++) thrds[i].thrd.join(0);
			double sum = 0;
			rowSums = new double[n+1];
			colSums = new double[n+1];
			for (int i=0; i<n; i++) {
				for (int k=0; k<n; k++) {
					sum += scores[i][k];
					rowSums[i] += scores[i][k];
					colSums[k] += scores[i][k];
				}
			}
			for (int i=0; i<n; i++) {
				scores[n][i] = sum/colSums[i];
				scores[i][n] = sum/rowSums[i];
				rowSums[n] += scores[n][i];
				colSums[n] += scores[i][n];
			}
			
			for (int i=0; i<n+1; i++) {
				for (int k=0; k<n+1; k++) {
					in_scores[i][k] = scores[i][k] / (colSums[k]+scores[n][k]);
					out_scores[i][k] = scores[i][k] / (rowSums[i]+scores[i][n]);
					if (i<n && k<n) scores2[i][k] = scores[i][k];
				}
			}

			workloads = new int[threadNr][][];
			for (int i=0; i<threadNr; i++) {
				int low=n*n*i/threadNr, up=n*n*(i+1)/threadNr;
				workloads[i] = new int[up-low][2];
				for (int k=0; k<workloads[i].length; k++) {
					workloads[i][k][0] = (low+k)/n;
					workloads[i][k][1] = (low+k)%n;
				}
				thrds[i] = new ScoringQueue(workloads[i],false,true,false,i+1);
				thrds[i].start();
			}
			for (int i=0; i<threadNr; i++) thrds[i].thrd.join(0);
			
			double[][] tmp = new double[n][n];
			for (int i=0; i<n; i++) {
				for (int k=0; k<n; k++) {
					tmp[i][k] = scores2[i][k] * scores3[i][k];
					tmp[i][k] = Math.atan(tmp[i][k]);
					tmp[i][k] /= Math.PI/2;
					tmp[i][k] = Math.pow(n*n,tmp[i][k]);
					scores[i][k] *= tmp[i][k];
				}
			}
			
			sum = 0;
			rowSums = new double[n+1];
			colSums = new double[n+1];
			for (int i=0; i<n; i++) {
				for (int k=0; k<n; k++) {
					sum += scores[i][k];
					rowSums[i] += scores[i][k];
					colSums[k] += scores[i][k];
					scores2[i][k] = scores3[i][k];
				}
			}
			for (int i=0; i<n; i++) {
				scores[n][i] = sum/colSums[i];
				scores[i][n] = sum/rowSums[i];
			}
			
			for (int i=0; i<n+1; i++)
				for (int k=0; k<n+1; k++) {
					scores3[i][k] = scores[i][k];
					int sign = 0;
					if (scores3[i][k]>0) sign = 1;
					else if (scores3[i][k]<0) sign = -1;
					scores3[i][k] = Math.exp(sign * Math.log(Math.abs(scores3[i][k]+1)));
				}
			
			rowSums = new double[n+1];
			colSums = new double[n+1];
			
			for (int i=0; i<n+1; i++)
				for (int k=0; k<n+1; k++) {
					rowSums[i] += scores3[i][k];
					colSums[k] += scores3[i][k];
				}
			
			for (int i=0; i<n+1; i++) 
				for (int k=0; k<n+1; k++) 
					scores3[i][k] = scores3[i][k]*scores3[i][k] / rowSums[i] / colSums[k];
			
			ArrayList<int[]> mergepos = new ArrayList<int[]>();
			for (int i=0; i<n; i++)
				for (int k=0; k<n; k++) {
					int[] pos = {i,k};
					if (scores3[i][k] > confidence_threshold) {
						int j=0;
						for (; j<mergepos.size(); j++) {
							int[] pos2 = mergepos.get(j);
							if (scores3[pos2[0]][pos2[1]]<scores3[pos[0]][pos[1]]) break;
						}
						mergepos.add(j,pos);
					}
				}
			
			for (int i=mergepos.size()-1; i>-1; i--) 
				for (int k=0; k<mergepos.size(); k++) {
					if (mergepos.get(i)[0] == mergepos.get(k)[0] &&
							scores3[mergepos.get(i)[0]][mergepos.get(i)[1]]<scores3[mergepos.get(k)[0]][mergepos.get(k)[1]]) {
						mergepos.remove(i);
						break;
					}
					if (mergepos.get(i)[1] == mergepos.get(k)[1] &&
							scores3[mergepos.get(i)[0]][mergepos.get(i)[1]]<scores3[mergepos.get(k)[0]][mergepos.get(k)[1]]) {
						mergepos.remove(i);
						break;
					}
					if (mergepos.get(i)[0] == mergepos.get(k)[1] && mergepos.get(i)[1] == mergepos.get(k)[0] &&
							scores3[mergepos.get(i)[0]][mergepos.get(i)[1]]<scores3[mergepos.get(k)[0]][mergepos.get(k)[1]]) {
						mergepos.remove(i);
						break;
					}
				}
			
			if (mergepos.size()==0)
				if (force_result) {
					if (print_progbar) System.out.print("Not too sure about the next merge...\r\n");
					order = new double[n][n];
					progress = 0;
					percent = 0;
					
					for (int i=0; i<threadNr; i++) {
						thrds[i] = new ScoringQueue(workloads[i],false,false,true,i+1);
						thrds[i].start();
					}
					for (int i=0; i<threadNr; i++) thrds[i].thrd.join(0);
					
					double max = 0;
					int[] maxi = {0,0};
					for (int i=0; i<n; i++)
						for (int k=0; k<n; k++) 
							if (order[i][k] > max) {
								max = order[i][k];
								maxi[0] = i;
								maxi[1] = k;
							}
					
					mergepos.add(maxi);
				} else {
					System.out.println("Found no clear match... aborting assembly");
					abort = true;
					break;
				}
			
			for (int i=0; i<mergepos.size(); i++)
				for (int k=i+1; k<mergepos.size(); k++)
					if (mergepos.get(i)[1]>mergepos.get(k)[1]) {
						
					}
			
			ArrayList<Integer> delete_list = new ArrayList<Integer>();
			for (int i=0; i<mergepos.size(); i++) {
				int a = mergepos.get(i)[0];
				int b = mergepos.get(i)[1];
				//if (ignore_list.contains(a) || ignore_list.contains(b)) continue;
				//ignore_list.add(a);
				//ignore_list.add(b);
				for (int k=i+1; k<mergepos.size(); k++) {
					int c = mergepos.get(k)[0];
					int d = mergepos.get(k)[1];
					if (c==b) {
						mergepos.get(k)[0] = a;
						scores3[a][d] = scores3[c][d];
						scores[a][d] = scores[c][d];
					}
					if (d==b) {
						mergepos.get(k)[1] = a;
						scores3[c][a] = scores3[c][d];
						scores[c][a] = scores[c][d];
					}
					if (c==b && d==a) mergepos.remove(k);
				}

				reads.set(a,overlap(reads.get(a),reads.get(b)));
				boolean appended = false;
				for (int k=0; k<delete_list.size(); k++)
					if (delete_list.get(k)<b) {
						delete_list.add(k, b);
						appended = true;
						break;
					}
				if (!appended) delete_list.add(b);
			}
			
			for (int i : delete_list) 
				reads.remove(i);
			delete_substring_reads();
			n = reads.size();
			System.out.println(" finished with loop");
		}
		if (!abort) {
			String assembly;
			
			if (n==2) {
				scores = new double[n+1][n+1];
				merge_results = new String[n][n];
				score(0,1);
				int start = 0;
				if (scores[0][1]<scores[1][0]) start = 1;
				assembly = overlap(reads.get(start),reads.get(1-start));
			} else assembly = reads.get(0);
			
			System.out.println();
			System.out.printf("\n             Final assembly : \n%s",assembly);
			if (print_alignment) System.out.println("\nAlignment of true sequence and assembled sequence:");
			new ErrorCount(truth,assembly,line_width,print_alignment);
		} else {
			System.out.println("               Leftover reads:\n");
			for (String r : reads) System.out.println(r);
		}
		System.out.println();
		long t = System.currentTimeMillis()-t0;
		System.out.printf("                    Runtime : \t%02d:%02d:%02d\n\n",t/3600000,(t%3600000)/60000,(t%60000)/1000);
		
	}
	
	public static void get_order(int a, int b) {
		for (int r=0; r<n; r++)
			for (int c=0; c<n; c++) { 		// n^4 loop =( 				// weighing the different matrices differently
				if (scores2[r][c] < scores2[a][b]) order[a][b] += 1;	// weight of 'better-to-merge-this-first-matrix' 
				if (scores3[r][c] < scores3[a][b]) order[a][b] += 2.5; 	// weight of 'nothing-comes-between-these-two-matrix'
			}
	}
	
	public static void score(int a, int b) {
		
		String seq1=","+reads.get(a),
			   seq2=","+reads.get(b);
		int l1=seq1.length(),l2=seq2.length();
		double[][] sm = new double[l1][l2];
		int[][][]  pm = new int[l1][l2][2];
		int r1=l1-1,c1=l2-1,r2=l1-1,c2=l2-1;
		
		for (int i=1; i<l1; i++)
			for (int k=1; k<l2; k++) {
				double x = 0;
				if (seq1.charAt(i)!=seq2.charAt(k)) x=sub;
				else x=-1; // added
				if (sm[i-1][k-1]+x <= sm[i-1][k]+gap && sm[i-1][k-1]+x <= sm[i][k-1]+gap) {
					sm[i][k] = sm[i-1][k-1]+x;
					pm[i][k][0] = i-1;
					pm[i][k][1] = k-1;
				} else if (sm[i-1][k]+gap <= sm[i-1][k-1]+x && sm[i-1][k] <= sm[i][k-1]) {
					sm[i][k] = sm[i-1][k]+gap;
					pm[i][k][0] = i-1;
					pm[i][k][1] = k;
				} else if (sm[i][k-1]+gap <= sm[i-1][k-1]+x && sm[i][k-1] <= sm[i-1][k]) {
					sm[i][k] = sm[i][k-1]+gap;
					pm[i][k][0] = i;
					pm[i][k][1] = k-1;
				}
			}
		for (int k=1,i=l1-1; k<l1 && k<l2; k++) if (sm[i][k] <= sm[i][c1]) c1 = k;
		for (int i=1,k=l2-1; i<l2 && i<l1; i++) if (sm[i][k] <= sm[r2][k]) r2 = i;
		
		String ol11="",ol12="",
			   ol21="",ol22="";
		ArrayList<int[]> trace1 = new ArrayList<int[]>(),
						 trace2 = new ArrayList<int[]>();
		int[] abmin={r1,c1}, bamin={r2,c2};
		
		while (c1>0) {
			int[] p = pm[r1][c1];
			if (p[1]==c1) ol12 = "_" + ol12;
			else ol12 = seq2.charAt(c1) + ol12;
			if (p[0]==r1) ol11 = "_" + ol11;
			else ol11 = seq1.charAt(r1) + ol11;
			int[] tmp = {1,1};
			if (r1>1) tmp[0] = r1;
			if (c1>1) tmp[1] = c1;
			trace1.add(tmp);
			r1 = p[0];
			c1 = p[1];
		}
		for (int i=ol11.length()-1; i>-1; i--) {
			if (ol11.charAt(i)=='_' || ol12.charAt(i)=='_') {
				// insnr1++;
				ol11 = ol11.substring(0,i) + ol12.charAt(i) + ol11.substring(i+1);
			}
		}
		
		while (c2>0) {
			int[] p = pm[r2][c2];
			if (p[1]==c2) ol22 = "_" + ol22;
			else ol22 = seq2.charAt(c2) + ol22;
			if (p[0]==r2) ol21 = "_" + ol21;
			else ol21 = seq1.charAt(r2) + ol21;
			int[] tmp = {1,1};
			if (r2>1) tmp[0] = r2;
			if (c2>1) tmp[1] = c2;
			trace2.add(tmp);
			r2 = p[0];
			c2 = p[1];
		}
		for (int i=ol21.length()-1; i>-1; i--) {
			if (ol21.charAt(i)=='_' || ol22.charAt(i)=='_') {
				// insnr2++;
				ol21 = ol21.substring(0,i) + ol22.charAt(i) + ol21.substring(i+1);
			}
		}
		
		String newseq1 = seq1.substring(1,trace1.get(trace1.size()-1)[0]) + ol11 + seq2.substring(trace1.get(0)[1]+1);
		String newseq2 = seq2.substring(1,trace2.get(trace2.size()-1)[1]) + ol21 + seq1.substring(trace2.get(0)[0]+1);
		
		merge_results[a][b] = newseq1;
		merge_results[b][a] = newseq2;
		scores[a][b] = Math.pow(alphabet_size,-sm[abmin[0]][abmin[1]]);//Math.exp(-sm[abmin[0]][abmin[1]]) / (1.-Math.pow(1.-1./alphabet_size,insnr1+1)) * Math.pow(alphabet_size,trace1.size()-1-insnr1);
		scores[b][a] = Math.pow(alphabet_size,-sm[bamin[0]][bamin[1]]);//Math.exp(-sm[bamin[0]][bamin[1]]) / (1.-Math.pow(1.-1./alphabet_size,insnr2+1)) * Math.pow(alphabet_size,trace2.size()-1-insnr2);
	}
	
	public static void deep_score(int a,int b) {
		if (a==b) {
			return;
		}
		
		double outSum = 0, inSum = 0;
		
		for (int j=0;j<n;j++) {
			if (j==a || j==b) continue;
			String seq1=","+merge_results[a][b],
				   seq2=","+reads.get(j);
			
			int l1=seq1.length(),l2=seq2.length();
			double[][] sm = new double[l1][l2];
			int[][][]  pm = new int[l1][l2][2];
			int r1=l1-1,c1=l2-1,r2=l1-1,c2=l2-1;
			
			for (int i=1; i<l1; i++)
				for (int k=1; k<l2; k++) {
					double x = 0;
					if (seq1.charAt(i)!=seq2.charAt(k)) x=sub;
					else x=-1; // added
					if (sm[i-1][k-1]+x <= sm[i-1][k]+gap && sm[i-1][k-1]+x <= sm[i][k-1]+gap) {
						sm[i][k] = sm[i-1][k-1]+x;
						pm[i][k][0] = i-1;
						pm[i][k][1] = k-1;
					} else if (sm[i-1][k]+gap <= sm[i-1][k-1]+x && sm[i-1][k] <= sm[i][k-1]) {
						sm[i][k] = sm[i-1][k]+gap;
						pm[i][k][0] = i-1;
						pm[i][k][1] = k;
					} else if (sm[i][k-1]+gap <= sm[i-1][k-1]+x && sm[i][k-1] <= sm[i-1][k]) {
						sm[i][k] = sm[i][k-1]+gap;
						pm[i][k][0] = i;
						pm[i][k][1] = k-1;
					}
				}
			for (int k=1,i=l1-1; k<l1 && k<l2; k++) if (sm[i][k] <= sm[i][c1]) c1 = k;
			for (int i=1,k=l2-1; i<l2 && i<l1; i++) if (sm[i][k] <= sm[r2][k]) r2 = i;
			
			double score1 = Math.pow(4,-sm[r1][c1]);
			double score2 = Math.pow(4,-sm[r2][c2]);
			outSum += score1;
			inSum  += score2;
			score1 *= out_scores[a][j] * in_scores[a][j];
			score2 *= out_scores[j][b] * in_scores[j][b];
			
			// sequence of strings :: score1 : (a->b)->j  ;; score2 : j->(a->b)
			// -> if a->j is likely, then subtract from that likelihood all scores a->x->j
			// but all x->a->j combs are added to a->j, so if a->j is a good merge, it will come out on top
			
			// can i parallelize this so it can run on a GPU? (problem: volatile arrays...) 
			scores2[a][j] -= score1;
			scores2[j][b] -= score2;
			scores2[b][j] += score1;
			scores2[j][a] += score2;
		}
		
		outSum += scores[b][n];
		inSum += scores[n][a];
		
		double ratio = outSum / (rowSums[b]+scores[b][n]);
		ratio *= inSum / (colSums[a]+scores[n][a]);
		scores3[a][b] = ratio;
	}
	
	public static String overlap (String seq1,String seq2){
		seq1=","+seq1;
		seq2=","+seq2;
		int l1=seq1.length(),l2=seq2.length();
		double[][] sm = new double[l1][l2];
		int[][][]  pm = new int[l1][l2][2];
		int r=l1-1,c=l2-1;
		int[] end = {r,c};
		for (int i=1; i<l1; i++)
			for (int k=1; k<l2; k++) {
				double x = 0;
				if (seq1.charAt(i)!=seq2.charAt(k)) x=sub;
				else x=-1; // added
				if (sm[i-1][k-1]+x <= sm[i-1][k]+gap && sm[i-1][k-1]+x <= sm[i][k-1]+gap) {
					sm[i][k] = sm[i-1][k-1]+x;
					pm[i][k][0] = i-1;
					pm[i][k][1] = k-1;
				} else if (sm[i-1][k]+gap <= sm[i-1][k-1]+x && sm[i-1][k] <= sm[i][k-1]) {
					sm[i][k] = sm[i-1][k]+gap;
					pm[i][k][0] = i-1;
					pm[i][k][1] = k;
				} else if (sm[i][k-1]+gap <= sm[i-1][k-1]+x && sm[i][k-1] <= sm[i-1][k]) {
					sm[i][k] = sm[i][k-1]+gap;
					pm[i][k][0] = i;
					pm[i][k][1] = k-1;
				}
			}
		
		for (int k=1,i=l1-1; k<l2; k++) if (sm[i][k] <= sm[i][c]) c = k;
		for (int i=1,k=l2-1; i<l1; i++) if (sm[i][k] <= sm[r][k]) r = i;
		if (sm[l1-1][c]<sm[r][l2-1]) {
			end[0] = l1-1; end[1] = c;
		} else {
			end[0] = r; end[1] = l2-1;
		}
		r = end[0]; c = end[1];
		String ol1="",ol2="";
		ArrayList<int[]> trace = new ArrayList<int[]>();
		while (c>0 && r>0) {
			int[] p = pm[r][c];
			if (p[1]==c) ol2 = "_" + ol2;
			else ol2 = seq2.charAt(c) + ol2;
			if (p[0]==r) ol1 = "_" + ol1;
			else ol1 = seq1.charAt(r) + ol1;
			int[] tmp = {1,1};
			if (r>1) tmp[0] = r;
			if (c>1) tmp[1] = c;
			trace.add(tmp);
			r = p[0];
			c = p[1];
		}
		
		seq1 = spaces(c) + seq1.substring(1,r+1) + ol1 + seq1.substring(end[0]+1);
		seq2 = spaces(r) + seq2.substring(1,c+1) + ol2 + seq2.substring(end[1]+1);
		if (seq1.length()<seq2.length()) seq1 += spaces(seq2.length()-seq1.length());
		if (seq2.length()<seq1.length()) seq2 += spaces(seq1.length()-seq2.length());
		String al = spaces(Math.max(r, c)), res="";
		for (int i=0; i<seq1.length(); i++) {
			if (seq1.charAt(i) != ' '  && seq2.charAt(i) != ' ') {
				if (seq1.charAt(i)==seq2.charAt(i)) al += "|"; 
				else al += "*";
			}
			if (seq1.charAt(i) != ' ' || seq2.charAt(i) != ' ') {
				if (seq1.charAt(i) != '_' && seq1.charAt(i) != ' ') res += seq1.charAt(i);
				else res += seq2.charAt(i);
			}
		}
		al += spaces(seq1.length()-al.length());
		if (seq1.length() > line_width) {
			int alstart = r+c, alend = r+c+ol1.length();
			int alcenter = (alstart+alend)/2;
			boolean left1=false,leftal=false,left2=false,
					right1=false,rightal=false,right2=false;
			while (alcenter>line_width/2-3) {
				if (seq1.charAt(0)!=' ') left1=true;
				if (al.charAt(0)!=' ') leftal=true;
				if (seq2.charAt(0)!=' ') left2=true;
				seq1 = seq1.substring(1);
				al = al.substring(1);
				seq2 = seq2.substring(1);
				alcenter--;
				alstart--;
				alend--;
			}
			while (seq1.length()>line_width) {
				if (seq1.charAt(seq1.length()-1)!=' ') right1=true;
				if (al.charAt(al.length()-1)!=' ') rightal=true;
				if (seq2.charAt(seq2.length()-1)!=' ') right2=true;
				seq1 = seq1.substring(0,seq1.length()-1);
				al = al.substring(0,al.length()-1);
				seq2 = seq2.substring(0,seq2.length()-1);
			}
			
			if (left1) seq1 = "..." + seq1;
			else seq1 = "   " + seq1;
			if (right1) seq1 += "...";
			if (leftal) al = "..." + al;
			else al = "   " + al;
			if (rightal) al += "...";
			if (left2) seq2 = "..." + seq2;
			else seq2 = "   " + seq2;
			if (right2) seq2 += "...";
			
		}
		if (print_merges) System.out.printf("\n\t%s\n\t%s\n\t%s\n",seq1,al,seq2);
		return res;
	}
	
	public static String spaces(int n) {
		return new String(new char[n]).replace("\0"," ");
	}
	
	static class ScoringQueue implements Runnable {
		
		private boolean scoring_,deep_scoring_,counting_;
		private int[][] indices_;
		private Thread thrd;
		
		String threadID;
		
		ScoringQueue (int[][] indices, boolean scoring, boolean deep_scoring, boolean counting, int thread) {
			scoring_ = scoring;
			deep_scoring_ = deep_scoring;
			counting_ = counting;
			indices_ = indices;
			threadID = String.valueOf(thread);
		}

		@Override
		public void run() {
			if (scoring_) {
				for (int[] pos : indices_) {
					score(pos[0],pos[1]);
				}
			} else if (deep_scoring_) {
				for (int[] pos : indices_) {
					deep_score(pos[0],pos[1]);
					progress++;
					if (print_progbar && (int)(progress*100/(n*n))==percent) {
						System.out.print(":\r\n");
						percent++;
					}
				}
			} else if (counting_){
				for (int[] pos : indices_) {
					get_order(pos[0],pos[1]);
					progress++;
					if (print_progbar && (int)(progress*100/(n*n))==percent) {
						System.out.print(".\r\n");
						percent++;
					}
				}
			} else {
				int[] indices = indices_[0];
				char[] alph = new char[alphabet.size()+1];
				alph[0] = ',';
				for (int i=1; i<alph.length; i++)
					alph[i] = alphabet.get(i-1);
				for (int index : indices) {
					new ErrorCorrection(reads,index,alph,coverage_estimate,subchance,gapchance);
					progress++;
					if (print_progbar && (int)(progress*100/m)==percent) {
						System.out.print("~\r\n");
						percent++;
					}
				}
			}
			
		}
		
		public void start() {
			thrd = new Thread(this,threadID);
			thrd.start();
		}
	}

	
}
