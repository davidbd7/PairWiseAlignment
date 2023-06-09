import java.util.*;
import java.io.*;

public class PairwiseAlignment {
	
	public static void main(String[] args){
		//creates instance of class PWA
		PairwiseAlignment pwa = new PairwiseAlignment();
		
		//creates new file
		File f = new File("fastafile.txt");
		
		//read the file
		ArrayList<String> sequences = pwa.readFile(f);
		
		//run the alignment alg
		ArrayList<int[][]> matrixList = pwa.align(sequences.get(0), sequences.get(1));
		
		//print resulting matrices
		pwa.printMatrix(matrixList.get(0));
		pwa.printMatrix(matrixList.get(1));
		
		//print the score
		System.out.println("Score of the Alingment " + pwa.scoreOfMatrix(matrixList.get(0)) + "\n");
		
		//print the trace back
		String[] alignedStrings = pwa.traceback(matrixList.get(1), sequences.get(0), sequences.get(1));
		System.out.println("First Sequence: " + alignedStrings[0] + "\n");
		System.out.println("Second Sequence: " + alignedStrings[1] + "\n");
	}
	
	public void printMatrix(int[][] m) {
		//prints the given matrix
		
		//traverse the matrix
		for(int i = 0; i < m.length; i++) {
			for(int j = 0; j < m[0].length; j++) {
				//print out the value at i,j
				System.out.print(m[i][j] + "\t");
			}
			//print new line
			System.out.println("");
		}
		
		//print new line for format
		System.out.println();
	}
	
	public String arrayListToString(ArrayList<Character> l) {
		StringBuilder builder = new StringBuilder(l.size());
		for(int i = 0; i < l.size(); i++) {
			builder.append(l.get(i));
		}
		return builder.toString();
	}

	public ArrayList<String> readFile(File f){
		//create an arraylist to return the results
		ArrayList<String> result = new ArrayList<>();
		
		//make sure the file is in fasta format
		try{
			Scanner scan = new Scanner(f);
			while(scan.hasNextLine()) {
				//read the def line
				scan.nextLine();
				//read the sequence
				result.add(scan.nextLine());
			}
			
		}catch(FileNotFoundException e) {
			System.out.println(e.getMessage());
		}
		return result;
	}

	//initializes the pairwise matrix
	public int[][] initScoreMatrix(int[][] m){
		//initialize 0,0 to 0
		m[0][0] = 0;
		
		//initialize all rows on the first column to -2 of the previous row
		for(int i = 1; i < m.length; i++) {
			m[i][0] = m[i - 1][0] - 2;
		}
		
		//initialize all columns on the first row to -2 of the previous column
		for(int i = 1; i < m[0].length; i++) {
			m[0][i] = m[0][i - 1] - 2;
		}
		
		//return the initialized matrix
		return m;
	}
	
	public int[][] initDirMatrix(int[][] m) {
		//0 = points to no direction. 1 = north. 2 = west. 3 = northwest
		
		//initialize 0,0 to 0 (no direction)
		m[0][0] = 0;
		
		//initialize the first column to point north
		for(int i = 1; i < m.length; i++) {
			m[i][0] = 1;
		}
		
		//initialize the first row to point west
		for(int i = 1; i < m[0].length; i++) {
			m[0][i] = 2;
		}
		
		//return the initialized matrix
		return m;
	}
	
	public int[] calculateScore(char c1, char c2, int north, int west, int northwest){
		//stores the score at index 0 and the direction at index 1
		int[] results = new int[2];
		
		//subtract 2 from the value to the north
		north -= 2;
		
		//subtract 2 from the value to the west
		west -= 2;
		
		//for the value to NW, compare the two characters, if they are equal +1, if they are not equal - 1
		if(c1 == c2) {
			northwest++;
		} else {
			northwest--;
		}
		
		//compare the three values and take the highest one
		results[0] = Math.max(northwest, Math.max(north, west));
		
		//determine direction. if there are two directions, just take the first one
		if(results[0] == north) {
			results[1] = 1;
		} else if(results[0] == west) {
			results[1] = 2;
		} else {
			results[1] = 3;
		}
		
		return results;
	}
	
	public ArrayList<int[][]> align(String first, String second){
		//create an arraylist to store the results
		ArrayList<int[][]> result = new ArrayList<>();
		
		//create two matrices of first.length + 1 by second.length() + 1
		
		//the first matrix is for keeping the score of the alignment
		int[][] score = new int[first.length() + 1][second.length() + 1];
		//the second matrix is for the direction the cell points to (N, W, NW)
		int[][] dir = new int[first.length() + 1][second.length() + 1];
		
		//initialize the score matrix
		score = initScoreMatrix(score);
		
		//initialize the direction matrix
		dir = initDirMatrix(dir);
		
		//align the two sequences using the pairwise alignment algorithm
		//traverse the matrix
		for(int i = 1; i < score.length; i++){
			for(int j = 1; j < score[i].length; j++){
				//call calulateScore for every position
				int[] resultOfScore = calculateScore(first.charAt(i - 1), second.charAt(j - 1), 
						score[i - 1][j], score[i][j - 1], score[i - 1][j - 1]);
				
				//set score at the index
				score[i][j] = resultOfScore[0];
				
				//set direction at the index
				dir[i][j] = resultOfScore[1];
			}
		}
		//return the matrices inside of the arraylist
		result.add(score);
		result.add(dir);
		return result;
	}

	public int scoreOfMatrix(int[][] m) {
		//compute the score of the pairwise alignment
		return m[m.length - 1][m[0].length - 1];
	}
	
	public String[] traceback(int[][] m, String first, String second){
		//array of strings containing resulting sequences
		String[] result = new String[2];
		
		//keeps track of matrix position
		int i = m.length - 1;
		int j = m[0].length - 1;
		
		//create two char ArraysLists for the results
		ArrayList<Character> resultN = new ArrayList<>();
		ArrayList<Character> resultW = new ArrayList<>();
		
		//keeps track of char in the strings
		int bufferN = 1;
		int bufferW = 1;
		
		//keeps track of the number of iterations
		int n = 1;
		
		//compute one of the possible trace backs
		while(m[i][j] != 0) {
			//check which direction the index points to
			if(m[i][j] == 1) {
				//north, thus reduce a row
				i--;
				//also introduce an indel to the string to the north (second)
				resultN.add('-');
				resultW.add(first.charAt(first.length() - bufferN));
				
				//reduce buffer for N
				bufferN++;
			} else if(m[i][j] == 2) {
				//west, thus reduce a column
				j--;
				//also introduce an indel to the string to the west (first)
				resultN.add(second.charAt(second.length() - bufferW));
				resultW.add('-');
				
				//reduce buffer for W
				bufferW++;
			}else {
				//northwest, thus reduce both a row and a column
				i--;
				j--;
				//no indels added
				resultN.add(second.charAt(second.length() - bufferW));
				resultW.add(first.charAt(first.length() - bufferN));
				
				//increase the buffer
				bufferN++;
				bufferW++;
			}
			n++;
		}
		//reverse the ArrayList
		Collections.reverse(resultN);
		Collections.reverse(resultW);
		
		//convert the ArrayList to a string
		result[0] = arrayListToString(resultN);
		result[1] = arrayListToString(resultW);
		
		return result;
	}
}