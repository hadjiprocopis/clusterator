package ahp.org.Clusterator;

import java.lang.StringBuilder;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.Vector;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.Instance;
import net.sf.javaml.core.DenseInstance;

import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.model.Model;
import de.lmu.ifi.dbs.elki.data.model.OPTICSModel;
import de.lmu.ifi.dbs.elki.data.model.KMeansModel;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRange;
import de.lmu.ifi.dbs.elki.database.ids.DBIDIter;
import de.lmu.ifi.dbs.elki.database.ids.DBIDRef;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.database.relation.Relation;

import ahp.org.Utils.*;

public class ClusteringResult {
	private int	myNumClusters,
			myNumDims;
	private double	myCentroids[/*cluster*/][/*dim*/],
			myBorders[/*cluster*/][/*0:from,1:to*/][/*dim*/];
	/* prob of data i belonging to cluster j based on distance.
	  last entry is where it belongs (largest prob) and is an int
	  last-but-one entry is the cluster whose centroid is nearest and is an int */
	private	double	myMembership[/*cluster+2*/][/*dataindex*/] = null,
			myNumMembers[/*cluster*/] = null;
	// metrics identified by a name are one for each cluster and the last in the array is the overall
	private HashMap<String,double[]> myMetrics = null;

/*	public void recalculate(boolean force){
		if( (this.myMetrics==null) || (force==true) ){
			this.calculate_membership(this.xxx);
		}
	}
*/
	// constructor : from ELKI clustering  - fast optics
	public ClusteringResult(
		Clustering<OPTICSModel> c,
		Relation<NumberVector> rel,
		DBIDRange ids,
		// ignore these shits
		boolean shit1, boolean shit2
	){
		List<Cluster<OPTICSModel>> allClusters = c.getAllClusters();
		this.myNumClusters = allClusters.size();

		Cluster<OPTICSModel> clu = allClusters.get(0);
		NumberVector apointcoords;
		for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
			// To get the vector use:
			this.myNumDims = rel.get(it).getDimensionality();
			break;
		}
		this.myCentroids = new double[this.myNumClusters][this.myNumDims];
		this.myBorders = new double[this.myNumClusters][2][this.myNumDims];

		double left_border[][] = new double[this.myNumClusters][this.myNumDims];
		double right_border[][] = new double[this.myNumClusters][this.myNumDims];
		
		double min[] = new double[this.myNumDims],
			max[] = new double[this.myNumDims];
		double centroids_i[], left_border_i[], right_border_i[];
		double aval, mean[];
		int i, j, k, num_members_i;
		for(i=this.myNumClusters;i-->0;){
			clu = allClusters.get(i);
			centroids_i = this.myCentroids[i];
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
				apointcoords = rel.get(it);
				for(j=this.myNumDims;j-->0;){
					aval = apointcoords.doubleValue(j);
					left_border_i[j] =
					right_border_i[j] =
					centroids_i[j] = aval;
				}
				break;
			}
			for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
				apointcoords = rel.get(it);
				for(j=this.myNumDims;j-->0;){
					aval = apointcoords.doubleValue(j);
					centroids_i[j] = aval;
					if( aval < left_border_i[j] ){ left_border_i[j] = aval; }
					if( aval > right_border_i[j] ){ right_border_i[j] = aval; }
				}
			}
		}
		for(i=this.myNumClusters;i-->0;){
			centroids_i = this.myCentroids[i];
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			num_members_i = allClusters.get(i).size();
			for(j=this.myNumDims;j-->0;){
				this.myBorders[i][0][j] = left_border_i[j];
				this.myBorders[i][1][j] = right_border_i[j];
				centroids_i[j] /= num_members_i;
			}
		}
		this.sort_clusters_wrt_centroid();
	}
	// constructor : from ELKI clustering (Model, this is produced by DBSCAN)
	public ClusteringResult(
		Clustering<Model> c,
		Relation<NumberVector> rel,
		DBIDRange ids,
		// ignore this, <Model> and <KMeansModel> can not be distinguished by java which thinks the constructors with kmeans and dbscan are the same
		boolean is_elki_CRAP
	){
		List<Cluster<Model>> allClusters = c.getAllClusters();
		this.myNumClusters = allClusters.size();

		Cluster<Model> clu = allClusters.get(0);
		NumberVector apointcoords;
		for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
			// To get the vector use:
			this.myNumDims = rel.get(it).getDimensionality();
			break;
		}
		this.myCentroids = new double[this.myNumClusters][this.myNumDims];
		this.myBorders = new double[this.myNumClusters][2][this.myNumDims];

		double left_border[][] = new double[this.myNumClusters][this.myNumDims];
		double right_border[][] = new double[this.myNumClusters][this.myNumDims];
		
		double min[] = new double[this.myNumDims],
			max[] = new double[this.myNumDims];
		double centroids_i[], left_border_i[], right_border_i[];
		double aval, mean[];
		int i, j, k, num_members_i;
		for(i=this.myNumClusters;i-->0;){
			clu = allClusters.get(i);
			centroids_i = this.myCentroids[i];
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
				apointcoords = rel.get(it);
				for(j=this.myNumDims;j-->0;){
					aval = apointcoords.doubleValue(j);
					left_border_i[j] =
					right_border_i[j] =
					centroids_i[j] = aval;
				}
				break;
			}
			for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
				apointcoords = rel.get(it);
				for(j=this.myNumDims;j-->0;){
					aval = apointcoords.doubleValue(j);
					centroids_i[j] = aval;
					if( aval < left_border_i[j] ){ left_border_i[j] = aval; }
					if( aval > right_border_i[j] ){ right_border_i[j] = aval; }
				}
			}
		}
		for(i=this.myNumClusters;i-->0;){
			centroids_i = this.myCentroids[i];
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			num_members_i = allClusters.get(i).size();
			for(j=this.myNumDims;j-->0;){
				this.myBorders[i][0][j] = left_border_i[j];
				this.myBorders[i][1][j] = right_border_i[j];
				centroids_i[j] /= num_members_i;
			}
		}
		this.sort_clusters_wrt_centroid();
	}
	// constructor : from ELKI clustering (Kmeans Model - yes we need one for each shit)
	public ClusteringResult(
		Clustering<KMeansModel> c,
		Relation<NumberVector> rel,
		DBIDRange ids
	){
		List<Cluster<KMeansModel>> allClusters = c.getAllClusters();
		this.myNumClusters = allClusters.size();

		Cluster<KMeansModel> clu = allClusters.get(0);
		NumberVector apointcoords;
		for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
			// To get the vector use:
			this.myNumDims = rel.get(it).getDimensionality();
			break;
		}
		this.myCentroids = new double[this.myNumClusters][this.myNumDims];
		this.myBorders = new double[this.myNumClusters][2][this.myNumDims];

		double left_border[][] = new double[this.myNumClusters][this.myNumDims];
		double right_border[][] = new double[this.myNumClusters][this.myNumDims];
		
		double min[] = new double[this.myNumDims],
			max[] = new double[this.myNumDims];
		double centroids_i[], left_border_i[], right_border_i[];
		double aval, mean[];
		int i, j, k, num_members;
		for(i=this.myNumClusters;i-->0;){
			clu = allClusters.get(i);
			mean = clu.getModel().getMean();
			centroids_i = this.myCentroids[i];
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			for(j=this.myNumDims;j-->0;){
				centroids_i[j] = mean[j];
			}
			for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
				apointcoords = rel.get(it);
				for(j=this.myNumDims;j-->0;){
					aval = apointcoords.doubleValue(j);
					left_border_i[j] =
					right_border_i[j] = aval;
				}
				break;
			}
			for(DBIDIter it = clu.getIDs().iter(); it.valid(); it.advance()) {
				apointcoords = rel.get(it);
				for(j=this.myNumDims;j-->0;){
					aval = apointcoords.doubleValue(j);
					centroids_i[j] = aval;
					if( aval < left_border_i[j] ){ left_border_i[j] = aval; }
					if( aval > right_border_i[j] ){ right_border_i[j] = aval; }
				}
			}
		}
		for(i=this.myNumClusters;i-->0;){
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			for(j=this.myNumDims;j-->0;){
				this.myBorders[i][0][j] = left_border_i[j];
				this.myBorders[i][1][j] = right_border_i[j];
			}
		}
		this.sort_clusters_wrt_centroid();
	}
	// constructor: from a JavaML clustering
	public ClusteringResult(Dataset result[]){
		this.myNumClusters = result.length;
		this.myNumDims = result[0].instance(0).noAttributes();
		this.myCentroids = new double[this.myNumClusters][this.myNumDims];
		this.myBorders = new double[this.myNumClusters][2][this.myNumDims];

		double left_border[][] = new double[this.myNumClusters][this.myNumDims];
		double right_border[][] = new double[this.myNumClusters][this.myNumDims];
		
		double min[] = new double[this.myNumDims],
			max[] = new double[this.myNumDims];
		double centroids_i[], left_border_i[], right_border_i[];
		double aval;
		int i, j, k, num_members;
		Instance aninstance;
		for(i=this.myNumClusters;i-->0;){
			Dataset acluster = result[i];
			num_members = acluster.size();
			aninstance = acluster.instance(0);
			centroids_i = this.myCentroids[i];
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			for(j=this.myNumDims;j-->0;){
				left_border_i[j] =
				right_border_i[j] = aninstance.value(j);
				centroids_i[j] = 0.0;
			}
			for(k=num_members;k-->1;){
				aninstance = acluster.instance(k);
				for(j=this.myNumDims;j-->0;){
					aval = aninstance.value(j);
					centroids_i[j] = aval;
					if( aval < left_border_i[j] ){ left_border_i[j] = aval; }
					if( aval > right_border_i[j] ){ right_border_i[j] = aval; }
				}
			}
			for(j=this.myNumDims;j-->0;){
				centroids_i[j] /= num_members;
			}
		}
		for(i=this.myNumClusters;i-->0;){
			left_border_i = left_border[i];
			right_border_i = right_border[i];
			for(j=this.myNumDims;j-->0;){
				this.myBorders[i][0][j] = left_border_i[j];
				this.myBorders[i][1][j] = right_border_i[j];
			}
		}
		this.sort_clusters_wrt_centroid();
	}
	// constructor: from anything which gives us a breaks[cluster][dim][0/from,1/to]
	// breaks is like borders
	// for Ndim data : data[][]
	public ClusteringResult(
		double dat[/*dataindex*/][/*dim*/],
		double abreaks[/*acluster*/][/*from=0,to=1*/][/*adim*/]
	){
		this.myNumClusters = abreaks.length;
		this.myNumDims = dat[0].length;
		this.myCentroids = new double[this.myNumClusters][this.myNumDims];
		this.myBorders = new double[this.myNumClusters][2][this.myNumDims];

		int num_data = dat.length;
		
		double centroids_i[], left_border_i[], right_border_i[], aval[];
		int i, j, k, num_members;
		Instance aninstance;
		for(i=this.myNumClusters;i-->0;){
			for(j=this.myNumDims;j-->0;){
				this.myBorders[i][0][j] = abreaks[i][0][j];
				this.myBorders[i][1][j] = abreaks[i][1][j];
			}
		}

		for(i=this.myNumClusters;i-->0;){
			centroids_i = this.myCentroids[i];
			left_border_i = this.myBorders[i][0];
			right_border_i = this.myBorders[i][1];
			num_members = 0;
			for(k=num_data;k-->0;){
				aval = dat[k];
				if( Utils.inRange(aval, left_border_i, right_border_i) ){
					for(j=this.myNumDims;j-->0;){
						centroids_i[j] += aval[j];
					}
					num_members++;
					break; // can not belong to anybody else
				}
			}
			if( num_members > 0 ){
				for(j=this.myNumDims;j-->0;){ centroids_i[j] /= num_members; }
			}
		}
		this.sort_clusters_wrt_centroid();
	}
	// for 1D data
	public ClusteringResult(
		double dat[/*dataindex*/],
		double abreaks[/*acluster*/][/*from=0,to=1*/]
	){
		this.myNumClusters = abreaks.length;
		this.myNumDims = 1;
		this.myCentroids = new double[this.myNumClusters][this.myNumDims];
		this.myBorders = new double[this.myNumClusters][2][this.myNumDims];

		int num_data = dat.length;
		
		double centroids_i[], left_border_i, right_border_i, aval;
		int i, j, k, num_members;
		Instance aninstance;
		for(i=this.myNumClusters;i-->0;){
			this.myBorders[i][0][0] = abreaks[i][0];
			this.myBorders[i][1][0] = abreaks[i][1];
		}

		for(i=this.myNumClusters;i-->0;){
			centroids_i = this.myCentroids[i];
			left_border_i = this.myBorders[i][0][0];
			right_border_i = this.myBorders[i][1][0];
			num_members = 0;
			for(k=num_data;k-->0;){
				aval = dat[k];
				if( Utils.inRange(aval, left_border_i, right_border_i) ){
					centroids_i[0] += aval;
					num_members++;
					break; // can not belong to anybody else
				}
			}
			if( num_members > 0 ){
				centroids_i[0] /= num_members;
			}
		}
		this.sort_clusters_wrt_centroid();
	}
	// calcs prob of data item belonging to one of the clusters as the distance from centroid
	// Ndim case:
	private void calculate_membership(
		double dat[/*index*/][/*dim*/]
	){
		if( this.myMembership != null ){ System.err.println("ClusteringResult.java : calculate_membership() : has already been calculated. Will not proceed."); return; }
		int	i, j, num_clusters1 = this.myNumClusters+1,
			num_data = dat.length,
			num_dims = dat[0].length;

		this.myMembership = new double[this.myNumClusters+2][num_data];
		this.myNumMembers = new double[this.myNumClusters];

		double	abreak_min[], abreak_max[], aval[], amem, centroid_i[],
			membership_i[];
		int num_members;

		double dummy[] = new double[num_data];
		for(j=num_data;j-->0;){ dummy[j] = 0.0; } // i know
		this.myNumMembers = new double[this.myNumClusters];

		//myBorders[/*cluster*/][/*0:from,1:to*/][/*dim*/];
		//myCentroids[/*cluster*/][/*dim*/]
		// where the nearest cluster wrt centroid distance is:
		double	membership_last_but_one[] = this.myMembership[this.myNumClusters];
		// the cluster we belong based on borders:
		double	membership_last[] = this.myMembership[num_clusters1];
		for(i=this.myNumClusters;i-->0;){
			abreak_min = this.myBorders[i][0];
			abreak_max = this.myBorders[i][1];
			centroid_i = this.myCentroids[i];
			membership_i = this.myMembership[i];
			for(j=num_data;j-->0;){
				aval = dat[j];
				membership_i[j] = amem = Utils.euclidean_distance(centroid_i, aval);
				if( amem < dummy[j] ){
					dummy[j] = amem;
					membership_last_but_one[j] = i;
				}
			}
		}
		for(i=this.myNumClusters;i-->0;){
			abreak_min = this.myBorders[i][0];
			abreak_max = this.myBorders[i][1];
			num_members = 0;
			for(j=num_data;j-->0;){
				aval = dat[j];
				if( Utils.inRange(aval, abreak_min, abreak_max) ){
					// this is the cluster we are within its borders
					membership_last[j] = i;
					num_members++;
				}
			}
			this.myNumMembers[i] = num_members;
		}
	}
	// calcs prob of data item belonging to one of the clusters as the distance from centroid
	// 1D case:
	private void calculate_membership(
		double dat[/*index*/]
	){
		if( this.myMembership != null ){ System.err.println("ClusteringResult.java : calculate_membership() : has already been calculated. Will not proceed."); return; }
		int	i, j, num_clusters1 = this.myNumClusters+1,
			num_data = dat.length,
			num_dims = 1;

		this.myMembership = new double[this.myNumClusters+2][num_data];
		this.myNumMembers = new double[this.myNumClusters];

		double	abreak_min, abreak_max, aval, amem, centroid_i,
			membership_i[];
		int num_members;

		double dummy[] = new double[num_data];
		for(j=num_data;j-->0;){ dummy[j] = 0.0; } // i know
		this.myNumMembers = new double[this.myNumClusters];

		//myBorders[/*cluster*/][/*0:from,1:to*/][/*dim*/];
		//myCentroids[/*cluster*/][/*dim*/]
		// where the nearest cluster wrt centroid distance is:
		double	membership_last_but_one[] = this.myMembership[this.myNumClusters];
		// the cluster we belong based on borders:
		double	membership_last[] = this.myMembership[num_clusters1];
		for(i=this.myNumClusters;i-->0;){
			abreak_min = this.myBorders[i][0][0];
			abreak_max = this.myBorders[i][1][0];
			centroid_i = this.myCentroids[i][0];
			membership_i = this.myMembership[i];
			for(j=num_data;j-->0;){
				aval = dat[j];
				membership_i[j] = amem = Math.abs(centroid_i - aval);
				if( amem < dummy[j] ){
					dummy[j] = amem;
					membership_last_but_one[j] = i;
				}
			}
		}
		for(i=this.myNumClusters;i-->0;){
			abreak_min = this.myBorders[i][0][0];
			abreak_max = this.myBorders[i][1][0];
			num_members = 0;
			for(j=num_data;j-->0;){
				aval = dat[j];
				if( Utils.inRange(aval, abreak_min, abreak_max) ){
					// this is the cluster we are within its borders
					membership_last[j] = i;
					num_members++;
				}
			}
			this.myNumMembers[i] = num_members;
		}
	}
	public	HashMap<String,double[]> metrics(
		double dat[/*index*/]
	){
		if( this.myMetrics != null ){
			System.out.println("ClusteringResult.java : metrics() : calculations already done!");
			return this.myMetrics;
		}
		this.calculate_membership(dat);

		this.myMetrics = new HashMap<String,double[]>();

		int	i, j,
			num_data = dat.length,
			num_dims = 1;

		// density
		double density[] = new double[this.myNumClusters+1];
		this.myMetrics.put("density", density);
		double tot = 0.0;
		for(i=this.myNumClusters;i-->0;){
			density[i] = this.myNumMembers[i] / Utils.euclidean_distance(this.myBorders[i][1], this.myBorders[i][0]);
			tot += density[i];
		}
		density[this.myNumClusters] = tot / this.myNumClusters;
		return this.myMetrics;
	}
	private void sort_clusters_wrt_centroid(){
		int i, newindex;
		_CentroidPseudo cp[] = new _CentroidPseudo[this.myNumClusters];
		for(i=this.myNumClusters;i-->0;){
			cp[i] = new _CentroidPseudo(i,this.myCentroids[i]);
		}
		Arrays.sort(cp);
		double cloned_centroids[][] = this.myCentroids.clone();
		double cloned_borders[][][] = this.myBorders.clone();
		for(i=this.myNumClusters;i-->0;){
			newindex = cp[i].index;
			if( newindex != i ){
				this.myCentroids[i] = cloned_centroids[newindex];
				this.myBorders[i] = cloned_borders[newindex];
			}
		}
	}
	// it returns the breaks (the borders) as:
	//		myBorders[/*cluster*/][/*0:from,1:to*/][/*dim*/];
	// these is exactly as borders
	public	double[][][] breaks(){ return this.myBorders; }
	// returns the breaks assuming 1D as : [/*cluster*/][/*0:from,1:to*/]
	public double[][] breaks1d(){
		int L = this.myNumClusters, i;
		double ret[][] = new double[L][2];
		for(i=L;i-->0;){
			ret[i][0] = this.myBorders[i][0][0];
			ret[i][1] = this.myBorders[i][1][0];
		}
		return ret;
	}

	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("number of clusters: "+this.myNumClusters+" (in "+this.myNumDims+" dimensions):\n");
		int i, j;
		for(i=0;i<this.myNumClusters;i++){
			sb.append("cluster "+(i+1)+" has "+(this.myNumMembers==null?"N/A":this.myNumMembers[i])+" members:\n");
			sb.append("\tcentroid: [");
			for(j=0;j<this.myNumDims;j++){
				sb.append(this.myCentroids[i][j]+",");
			}
			sb.deleteCharAt(sb.length()-1);
			sb.append("]\n");
			sb.append("\tborders: [");
			//myBorders[/*cluster*/][/*0:from,1:to*/][/*dim*/];

			for(j=0;j<this.myNumDims;j++){
				sb.append("\n\t\t"
					+this.myBorders[i][0][j]
					+" to "
					+this.myBorders[i][1][j]
				);
			}
			sb.append("\n\t]\n");
		}
		if( this.myNumDims == 1 ){
			String Rc = "c(";
			for(i=0;i<this.myNumClusters;i++){
				Rc += this.myBorders[i][0][0] + ",";
			}
			Rc += this.myBorders[i-1][1][0]+")";
			sb.append("as R-vector : "+Rc);
		}
		return sb.toString();
	}
	public double[][] centroids(){ return this.myCentroids; }
	public double[][][] borders(){ return this.myBorders; }
	public void centroids(double[][] c){ this.myCentroids = c; }
	public void borders(double[][][] p){ this.myBorders = p; }
}
class _CentroidPseudo implements Comparable<_CentroidPseudo> {
	int index;
	double values[];
	public _CentroidPseudo(int ind, double c[]){ this.index = ind; this.values = c;  }
	public int compareTo(_CentroidPseudo another){
		double sum1 = Utils.euclidean_distance(values);
		double sum2 = Utils.euclidean_distance(another.values);
		return Double.compare(sum1, sum2);
	}
	public String toString(){
		return new String(this.index+" : "+Arrays.toString(this.values));
	}
}
