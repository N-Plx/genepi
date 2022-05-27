#include "cpp.h"
#include "hepevt.h"
#include "track_vars.h"
#include "dump_file.h"
#include"vertex.h"

void dump_file(int mode, int proc, double xsec)
{
  if(mode == 0)
    { 
      for(int i=2; i<hepevt.NHEP; i++)
	{
	  if(i!=3){
	    fprintf(ptr,"%5d %5d %5d %15.8f %15.8f %15.8f %15.8f\n",
		    hepevt.IDHEP[i], hepevt.JDAHEP[i][0], 
		    hepevt.JDAHEP[i][1], hepevt.PHEP[i][0], hepevt.PHEP[i][1], 
		    hepevt.PHEP[i][2],hepevt.PHEP[i][4]);
	  }
	}
    }
  else if(mode == 1)
    {
      if(proc ==0)
	{
	  fprintf(ptr,"%6d %6d %8.4f %8.4f %6d %6d\n",
		  trk.Ntracks-3, trk.TarA, trk.Q2, trk.nu, trk.Theli, trk.Bheli);
	  for(int i=2; i<trk.Ntracks; i++)
	    {
	      if(i!=3){      
		if(i==2){
		  fprintf(ptr,"%5d %5d %5d %5d %5d\n",
			  1, trk.Type[i], trk.Type[i], i-1, i-1);
		  fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n",  
			  trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i], sqrt(trk.E[i]*trk.E[i]-trk.P[i]*trk.P[i]));
		  fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n", 
			  0.,0.,-68.18,0.,0.);
		}
		if(i!=2){
		  fprintf(ptr,"%5d %5d %5d %5d %5d\n",
			  1, trk.Type[i], trk.Type[i], i-2, i-2);
		  fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n",  
			  trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i], sqrt(trk.E[i]*trk.E[i]-trk.P[i]*trk.P[i]));
		  fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n", 
			  0.,0.,-68.18,0.,0.);
		}
	      }
	    }
	}
      else if(proc == 1)
	{
	  fprintf(ptr,"%6d %6d %8.4f %8.4f %6d %6d\n",
		  trk.Ntracks-4, trk.TarA, trk.Q2, trk.nu, trk.Theli, trk.Bheli);
	  for(int i=2; i<trk.Ntracks; i++)
	    {
	      if(i!=3 && i!=6){      
		fprintf(ptr,"%5d %5d %5d %5d %5d\n",
			1, trk.Type[i], trk.Type[i], 1, 1);
		if(trk.Type[i]!=22 && trk.Type[i]!=11 && trk.Type[i]!=-11)fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n",  
										  trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i], sqrt(trk.E[i]*trk.E[i]-trk.P[i]*trk.P[i]));
		if(trk.Type[i]==22 && trk.E[i]*trk.E[i]-trk.P[i]*trk.P[i]>=0)fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n",  
										     trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i], sqrt(trk.E[i]*trk.E[i]-trk.P[i]*trk.P[i]));
		if(trk.Type[i]==22 && trk.E[i]*trk.E[i]-trk.P[i]*trk.P[i]<0)fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n",  
										    trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i], 0.);
		if(trk.Type[i]==-11 || trk.Type[i]==11)fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n",  
							       trk.Px[i], trk.Py[i], trk.Pz[i], trk.E[i], 0.00051);
		fprintf(ptr,"%15.8f %15.8f %15.8f %15.8f %15.8f\n", 
			0.,0.,-68.18,0.,0.);
	      }
	    }
	}
    }
  else if(mode == 2)
    {
      int index = 1;
      double vtx_x,vtx_y,vtx_z;
      vertex_xyz(&vtx_x, &vtx_y, &vtx_z);
      if(proc == 0)
	{
	  fprintf(ptr,"%d %d %d %d %d %d %lf %d %d %lf\n", hepevt.NHEP-3, 2, 1, 0, 0, 
		  11, hepevt.PHEP[0][3], hepevt.IDHEP[1], proc, xsec);
	  for(int i=2; i<hepevt.NHEP; i++)
	    {
	      if(i!=3){
		//Skip stupid events                                                                                                                                                                               
                if(hepevt.IDHEP[i] == -1000) continue;
		// index, lifetime, active status
		fprintf(ptr,"%d %d %d ", index, -1, 1); 
		index++;
		fprintf(ptr,"%5d %5d %5d %15.8f %15.8f %15.8f %15.8f %lf ",
			hepevt.IDHEP[i], 0, 0,
			hepevt.PHEP[i][0], hepevt.PHEP[i][1],
			hepevt.PHEP[i][2],hepevt.PHEP[i][3], hepevt.PHEP[i][4]);
		//vertex position
		fprintf(ptr,"%lf %lf %lf\n", vtx_x, vtx_y, vtx_z);
	      }
	    }
	}
      else if(proc == 1)
	{
	  fprintf(ptr,"%d %d %d %d %d %d %lf %d %d %lf\n", hepevt.NHEP-3-1, 2, 1, 0, 0, 
		  11, hepevt.PHEP[0][3], hepevt.IDHEP[1], proc, xsec);
	  for(int i=2; i<hepevt.NHEP; i++)
	    {
	      if(i!=3 && i!=6){
		// index, lifetime, active status
		fprintf(ptr,"%d %d %d ", index, -1, 1); 
		index++;
		fprintf(ptr,"%5d %5d %5d %15.8f %15.8f %15.8f %15.8f %lf ",
			hepevt.IDHEP[i], 0, 0,
			hepevt.PHEP[i][0], hepevt.PHEP[i][1],
			hepevt.PHEP[i][2],hepevt.PHEP[i][3], hepevt.PHEP[i][4]);
		// vertex: x, y, z position
		fprintf(ptr,"%lf %lf %lf\n", vtx_x, vtx_y, vtx_z);
	      }
	    }
	}      
    }
}
