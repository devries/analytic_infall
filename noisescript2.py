#!/usr/bin/env python
import sys
import os
import math
import numarray.random_array as rnd

def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv)!=6:
        print "Usage: %s <method> <frequency> <vmin> <vmax> <bestfitdatafile>"%argv[0]
        return 1

    hybrid = "%s_hybrid"%argv[1]
    simplex = "%s_simplex"%argv[1]
    dfile = argv[5]
    hybrid_out = "%s.hybrid_output"%dfile
    simplex_out = "%s.simplex_output"%dfile
    temp_in = "%s.temp_input"%dfile
    frequency = argv[2]
    vmin = argv[3]
    vmax = argv[4]

    fin = open(dfile,"r")
    params = []
    paramnames = []
    v = []
    tin = []
    tfit = []
    for l in fin:
        tokens = l.split()
        if tokens[0]=='#':
            if tokens[1]!="Attained":
                paramnames.append(tokens[1])
                params.append(float(tokens[2]))
        else:
            v.append(float(tokens[0]))
            tin.append(float(tokens[1]))
            tfit.append(float(tokens[2]))

    fin.close()
    
    tres = [tfit[i]-tin[i] for i in xrange(len(tin))]
    rms = reduce(lambda x,y: x+y*y, tres,0.0)   
    rms = math.sqrt(rms/len(tres))
    prmstr = ' '.join([str(i) for i in params])
    
    rnd.seed()
    meanparams = [0.0]*len(params)
    rmsparams = [0.0]*len(params)
    niter = 100
    for i in xrange(niter):
        ret = 1
        while(ret!=0):
            rmsarr = rnd.normal(0.0,rms,len(tfit))
            trand = [tfit[j]+rmsarr[j] for j in xrange(len(tfit))]
            fout = open(temp_in,"w")
            for j in range(len(v)):
                print >>fout,"%f\t%g"%(v[j],trand[j])
            fout.close()
            ret=os.system("%s %s %s %s %s %s %s"%(simplex,temp_in,frequency,vmin,vmax,prmstr,simplex_out))

        simplex_params = []
        fin = open(simplex_out,"r")
        for l in fin:
            tokens = l.split()
            if tokens[0]=='#':
                if tokens[1]!="Attained":
                    simplex_params.append(float(tokens[2]))
        fin.close()
        for j in xrange(len(params)):
            meanparams[j]+=simplex_params[j]
            rmsparams[j]+=simplex_params[j]*simplex_params[j]

    rmsparams=map(lambda x,y: math.sqrt((x-(y**2)/niter)/niter),rmsparams,meanparams)
    meanparams=[x/niter for x in meanparams]

    print "Residual RMS: ",rms
    print "Parameter\tMean\t\tSigma\t\tBest"
    for i in range(len(paramnames)):
        print "%s\t\t%f\t%f\t%f"%(paramnames[i],meanparams[i],rmsparams[i],params[i])

    os.remove(simplex_out)
    os.remove(temp_in)
    return 0
    
if __name__ == "__main__":
    sys.exit(main())
