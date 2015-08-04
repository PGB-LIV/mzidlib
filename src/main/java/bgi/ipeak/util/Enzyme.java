package bgi.ipeak.util;

public class Enzyme {	

	public final class Trypsin extends Enzyme {

		public Trypsin() {
			
		}
		public Boolean is_enzmyme_point(String n_termAA,String c_termAA) {
			return((n_termAA.equals("K")||n_termAA.equals("R"))&& !c_termAA.equals("P")||n_termAA.equals("-")||c_termAA.equals("-"));
		}
	}
	public enum EnzymeName{NO_ENZYME, TRYPSIN, CHYMOTRYPSIN, THERMOLYSIN, PROTEINASEK, PEPSIN, ELASTASE, 
	      LYSN, LYSC, ARGC, ASPN, GLUC
	      };
	EnzymeName enzymeName;
	static Enzyme theEnzyme=null;
	public Enzyme() {
		
	}
	
	
	
	public void SetEnzyme(EnzymeName enzyme_name) {
		if(enzyme_name ==EnzymeName.TRYPSIN){
			theEnzyme=new Trypsin();
		}
	}
}
