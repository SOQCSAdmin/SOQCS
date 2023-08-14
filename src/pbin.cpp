//======================================================================================================
// File pbin.cpp
//
// PROBABILITY BINS DEFINITION LIBRARY.
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "pbin.h"


//----------------------------------------
//
//  Creates a set of probability bins
//  The maximum number of photons and kets are set by default.
//
//----------------------------------------
p_bin::p_bin(int i_level):ket_list(i_level){
//  int i_level     // Number of levels to describe a ket in the set


    N=0;
    p=new double[maxket]();
}


//----------------------------------------
//
//  Creates a set of probability bins
//  The maximum number of kets is set by default.
//
//----------------------------------------
p_bin::p_bin(int i_nph, int i_level):ket_list(i_nph, i_level){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe a ket in the set


    N=0;
    p=new double[maxket]();
}


//-----------------------------------------------------
//
//  Creates a set of probability bins specifying the maximum number of bins.
//
//-----------------------------------------------------
p_bin::p_bin(int i_nph, int i_level, int i_maxket):ket_list(i_nph, i_level,i_maxket){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe a ket in the set
//  int i_maxket    // Maximum number of different kets in the set


    N=0;
    p=new double[maxket]();
}


//-----------------------------------------------------
//
//  Creates a set of probability bins specifying the maximum
//  number of bins  and the vector of equivalence between bin
//  levels and circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
p_bin::p_bin(int i_nph, int i_level, int i_maxket, int *i_vis):ket_list(i_nph, i_level,i_maxket, i_vis){
//  int i_nph       // Maximum number of photons
//  int i_level     // Number of levels to describe a ket in the set
//  int i_maxket    // Maximum number of different kets in the set
//  int *i_vis      // Vector of equivalence between state levels and circuit defined levels


    N=0;
    p=new double[maxket]();
}


//----------------------------------------
//
//  Destroys a state
//
//----------------------------------------
p_bin::~p_bin(){


    // Liberate memory of amplitude
    delete[] p;
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
p_bin *p_bin::clone(){
    // Variable
    p_bin *aux;     // Auxiliary state
    // Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index


    aux=new p_bin(nph, nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        aux->p[i]=p[i];
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }

    aux->ketindex=ketindex;
    aux->N=N;
    return aux;
}

//----------------------------------------
//
//  Clear a p_bin
//
//----------------------------------------
void p_bin::clear(){


    delete[] p;
    p=new double[maxket]();
    clear_kets();
}


//----------------------------------------
//
//  Counts a new sample
//
//----------------------------------------
int p_bin::add_count(int *occ){
//  int *occ;           // Occupation of those levels
//  Variables
    int index;          // Index to be returned


    index=add_ket(occ);
    if(index>=0) p[index]=p[index]+1.0;
    N=N+1;

    return index;
}


//----------------------------------------
//
//  Adds/sums the statistic from another
//  compatible bin set.
//
//----------------------------------------
int p_bin::add_bin(p_bin *input){
//  p_bin *input;     // Input probability bin set
//  Variables
    int    index;     // Index of each entry in this set
//  Auxiliary index
    int    i;         // Aux index


    for(i=0;i<input->nket;i++){
        index=add_ket(input->ket[i]);
        if(index<0) return -1;
        p[index]=p[index]+input->p[i];
    }

    N=N+input->N;
    return 0;
}


//----------------------------------------
//
// Adds/ sums the statistics form a quantum
// state computing the probability outcomes
// from it
//
//----------------------------------------
int p_bin::add_state(state *input){
//  state *input;     // Input state
//  Variables
    int    index;     // Index of each entry in this set
//  Auxiliary index
    int    i;         // Aux index


    for(i=0;i<input->nket;i++){
        index=add_ket(input->ket[i]);
        if(index<0) return -1;
        p[index]=p[index]+(double) abs(conj(input->ampl[i])*input->ampl[i]);
    }

    N=N+1;
    return 0;
}


//----------------------------------------
//
// Calculated the trace/total probability
// of the set of probability bins
//
//----------------------------------------
double p_bin::trace(){
//  Variables
    double M;      // Total probability
//  Auxiliary index
    int    i;      // Aux index


    M=0;
    for(i=0;i<nket;i++) M=M+p[i];
    return M/N;
}


//----------------------------------------
//
// Normalize a set of probability bins
// to one.
//
//----------------------------------------
void p_bin::normalize(){
//  Variables
    double M;      // Total probability
//  Auxiliary index
    int    i;      // Aux index

    M=0;
    for(i=0;i<nket;i++) M=M+p[i];
    for(i=0;i<nket;i++) p[i]=p[i]/M;
    N=1;
}


//----------------------------------------
//
// Returns the occupation of the indexed bin
// in string format.
//
//----------------------------------------
string p_bin::tag(int index){
//  int index;             // Index of the bin
//  Variables.
    int len;               // Length of the auxiliary string
    int diff;              // Difference with the number of levels
    long long int value;   // Hash value of the ket in that idex.
    string aux;            // Auxiliary string. Conversion of value
    string temp;           // Temporary string=aux+0's at the left to match nlevel length
//  Auxiliary index
    int i;                 // Aux index


    // Perform conversion
    value=decval(ket[index],nlevel,10);
    aux=to_string(value);
    // Match the nlevel lenfth with 0's at the left
    len=aux.length();
    diff=nlevel-len;
    for(i=0;i<diff;i++) temp.append("0",1);
    temp.append(aux);
    // Return result
    return temp;
}


//----------------------------------------
//
// Probability value of the indexed bin
//
//----------------------------------------
double p_bin::prob(int index){
//  int index;            // Index of the bin


    return p[index]/(double)N;
}


//----------------------------------------
//
// Probability value a bin from it definition
//
//----------------------------------------
double p_bin::prob(mati def,qocircuit *qoc){
//  mati       def;          // Definition of a probability bin in human readable form
//  qocircuit *qoc;          // Circuit where the detector is placed
//  Variables
    int index;               // Index in this set
    ket_list  *bra;          // Bra state created from the definition def


    // Create bra state
    bra= new ket_list(nph, nlevel,1,vis);
    bra->add_ket(def,qoc);
    index=find_ket(bra->ket[0]);
    delete bra;

    // Return probability value if present.
    if(index>=0) return prob(index);
    return 0.0;
}


//--------------------------------------------
//
// Prints the bin list
// ( various entries and probabilities )
//
//--------------------------------------------
void  p_bin::prnt_bins(){


    aux_prnt_bins(DEFFORMAT,0.0,false,nullptr);
}


//--------------------------------------------
//
// Prints the bin list
// ( various entries and probabilities )
//
//--------------------------------------------
void  p_bin::prnt_bins(double thresh){
//  double thresh;      // Probabilities below this value are not printed


    aux_prnt_bins(DEFFORMAT,thresh,false,nullptr);
}


//--------------------------------------------
//
// Prints the bin list
// ( in human readable form )
//
//--------------------------------------------
void  p_bin::prnt_bins(int format, double thresh, qocircuit *qoc){
//  int format;         // Notation of the bins to be printed
//  double thresh;      // Probabilities below this value are not printed
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)


    aux_prnt_bins(format, thresh,false,qoc);
}


//--------------------------------------------
//
// Prints the bin list
// ( in human readable form )
//
//--------------------------------------------
void p_bin::prnt_bins( int format, double thresh, bool loss, qocircuit *qoc){
//  int format;         // Notation of the bins to be printed
//  double thresh;      // Probabilities below this value are not printed
//  bool   loss;        // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the set of probability bins are referred. (Necessary to print in human readable form)


    aux_prnt_bins(format, thresh,loss,qoc);
}


//--------------------------------------------
//
// Prints the bin list
// ( in human readable form )
//  Uses qodev instead of qocircuit
//
//--------------------------------------------
void p_bin::prnt_bins(int format, double thresh, qodev *dev){
//  int format;         // Notation of the bins to be printed
//  double thresh;      // Probabilities below this value are not printed
//  qodev *dev;         // Device to which the set of probability bins are referred. (Necessary to print in human readable form)


    if(dev > (qodev *)nullptr) aux_prnt_bins(format, thresh, false, dev->circ);
    else aux_prnt_bins(format, thresh, false, nullptr);

}


//--------------------------------------------
//
// Prints the bin list
// ( in human readable form )
//  Uses qodev instead of qocircuit
//
//--------------------------------------------
void p_bin::prnt_bins( int format, double thresh, bool loss, qodev *dev){
//  int format;         // Notation of the bins to be printed
//  double thresh;      // Probabilities below this value are not printed
//  bool   loss;        // Print loss channels in different color? False=No/True=Yes
//  qodev *dev;         // Device to which the set of probability bins are referred. (Necessary to print in human readable form)


    if(dev > (qodev *)nullptr) aux_prnt_bins(format, thresh, loss ,dev->circ);
    else aux_prnt_bins(format, thresh, loss ,nullptr);
}


//--------------------------------------------
//
// Prints the bin list
// ( in human readable form )
//
//--------------------------------------------
void p_bin::aux_prnt_bins( int format, double thresh, bool loss, qocircuit *qoc){
//  int format;         // Notation of the bins to be printed
//  double thresh;      // Probabilities below this value are not printed
//  bool   loss;        // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the set of probability bins are referred. (Necessary to print in human readable form)
//  Variables
    int    firstline;   // Is this the first term? 1=Yes/0=No
    double prob;        // Probability
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        prob= p[i]/(double) N;
        if (prob>thresh){
            if(firstline==0){
                cout << endl;
            }
            firstline=0;

            cout <<right << setw(2) << i << " : ";
            prnt_ket(i,format,loss,qoc);
            cout << ": ";
            cout << left << setprecision(4) << prob;


        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//--------------------------------------------
//
//  Post selection process of the output state
//
//--------------------------------------------
p_bin *p_bin::post_selection(state *prj){
//  state *prj;       // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    int    index;     // Index of the sampled state
    int    npost;     // Number of levels in which post selection is performed
    int    selected;  // Is this state selected 0=No/1=Yes
    int   *islincl;   // Is level included? 0=No/1=Yes
    int   *occ;       // Occupation
    p_bin *nbin;      // New post-selected set of probabilities
// Auxiliary index
    int    i;         // Aux index
    int    j;         // Aux index
    int    k;         // Aux index
    int    l;         // Aux index


    // Determined the number of post-selected levels and which ones
    // will be included in nstate or not.
    npost=0;
    islincl=new int[nlevel];
    for(i=0;i<prj->nlevel;i++){
        if(prj->ket[0][i]>=0){
            islincl[i]=0;
            npost++;
        }else{
            islincl[i]=1;
        }
    }

    // Create new post-selected state / Reserve memory
    nbin=new p_bin(nph, nlevel-npost,maxket);
    k=0;
    for(l=0;l<nlevel;l++){
        if(islincl[l]==1){
            nbin->vis[k]=vis[l];
            k++;
        }
    }

    // For each projector term
    occ=new int[nlevel-npost]();
    for(i=0;i<prj->nket;i++){
        // Post-select each ket
        for(j=0;j<nket;j++){
            // Check selection condition
            selected=1;
            k=0;
            while((k<nlevel)&&(selected==1)){
                if ((ket[j][k]!=prj->ket[i][k])&&(prj->ket[i][k]>=0))  selected=0;
                k++;
            }

            // If is selected create the list of levels and
            // occupations for those not post-selected
            if(selected==1){
                k=0;
                for(l=0;l<nlevel;l++){
                    if(islincl[l]==1){
                        occ[k]=ket[j][l];
                        k++;
                    }
                }
                index=nbin->add_ket(occ);
                nbin->p[index]=nbin->p[index]+p[j]*pow(abs(prj->ampl[i]),2);
            }
        }
    }
    nbin->N=N;

    // Free memory
    delete[] occ;
    delete[] islincl;

    // Return bin
    return nbin;
}


//--------------------------------------------
//
//  Calculate a measurement from the stored statistics
//  using the corresponding circuit detector definitions
//
//--------------------------------------------
p_bin *p_bin::calc_measure(qocircuit *qoc){
//  qocircuit *qoc;         // Circuit with the detector definitions to calculate the measurement.
//  Variables
    int    S;               // Blinking and dark counts calculation number of iterations
    double stdev;           // Gaussian white noise standard deviation
    p_bin *inperiod;        // Output after removing photons that are not in the measurement period.
    p_bin *dark;            // Output after adding dark counts
    p_bin *blinked;         // Output after blinking effect considerations
    p_bin *lossed;          // Output after considering the detector efficiency
    p_bin *ignored;         // Output after removing ignored channels
    p_bin *measured;        // Output after considering detector conditions and removing loss channels (if needed)
    p_bin *counted;         // Output after removing temporal degrees of freedom (if needed)
    p_bin *noisy;           // Output after adding some gaussian white noise
    p_bin *aux;             // Auxiliary probability bin


    // Itialize variables
    S=qoc->R;
    stdev=sqrt(qoc->dev);

    // Compute Dark counts
    if(qoc->timed==0) dark=this->dark_counts(S,qoc);
    else dark=this->clone();

    // Compute Blinking
    blinked=dark->blink(S,qoc);

    // Measure
    if(qoc->losses==1) lossed=blinked->compute_loss(qoc);
    else lossed=blinked->clone();

    // Measurement window
    if(qoc->np>1) inperiod=lossed->meas_window(qoc);
    else inperiod=lossed->clone();

    // ignore channels
    if(qoc->nignored>0) ignored=inperiod->compute_ignored(qoc);
    else ignored=inperiod->clone();

    // Post-selection
    if(qoc->ncond>0) measured=ignored->compute_cond(qoc);
    else measured=ignored->clone();

    // Remove time
    if(qoc->ns>1){
       switch(qoc->timed){
            case 0:  // Counter
                aux=measured->perform_count(qoc);
                counted=aux->remove_time(qoc);
                delete aux;
                break;
            case 1:  // Clock: Partial trace
                counted=measured->remove_freq(qoc);
                break;
            case 2:  // Clock+Spectrum: Full
                counted=measured->clone();
                break;
            case 3:  // Clock: Manual mode
                counted=measured->remove_freq(qoc);
                break;
            case 4:  // Classify periods
                counted=measured->classify_period(qoc);
                break;
            default: // Default: Error
                cout << "calc_measure error: Time has a non-valid value" << endl;
                exit(0);
                break;

       }
    }else{
        counted=measured->clone();
    }

    // Add noise
    if (stdev>xcut) noisy=counted->white_noise(stdev);
    else noisy=counted->clone();

    // Free memory
    delete dark;
    delete blinked;
    delete lossed;
    delete inperiod;
    delete ignored;
    delete measured;
    delete counted;

    // Return final detection
    return noisy;
}

//--------------------------------------------
//
// Filters the photons out of the measurement window
// established by the detectors.
//
//--------------------------------------------
p_bin *p_bin::meas_window(qocircuit* qoc){
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    bool   filter;
    int    ch;          // Channel number
    int    m;           // Mode number
    int    s;           // Wavepacket number
    int    nwi;         // First period that can be measured
    int    nwf;         // Last period (+1) that can be measured
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // Set up measurement period.
    i=0;
    filter=false;
    while((filter==false)&&(i<qoc->nch)){
        if(qoc->det_win(0,i)>=0) filter=true;
        if(qoc->det_win(1,i)>=0) filter=true;
        i=i+1;
    }

    // Return bins if there is no filtering
    if(filter==false) return this->clone();

    // Calculate filtering
    newpbin=new p_bin(nph, nlevel,maxket,vis);
    for(i=0;i<nket;i++){
        occ=new int[nlevel]();
        for(ch=0;ch<qoc->ndetc;ch++){
            // Set filtering window for the current detector
            nwi=qoc->det_win(0,ch);
            nwf=qoc->det_win(1,ch)+1;
            if(qoc->det_win(0,ch)<0) nwi=0;
            if(qoc->det_win(1,ch)<0) nwf=qoc->np+1;

            // Store the occupation level if is inside the detection window
            for(m=0;m<qoc->nm;m++){
            for(s=0;s<qoc->ns;s++){
                l=qoc->i_idx[ch][m][s];
                // Store
                k=0;
                while(vis[k]!=l) k=k+1;
                if((s>=nwi*qoc->nsp)&&(s<nwf*qoc->nsp)) occ[k]=occ[k]+ket[i][k];

            }}
        }

        // Store occupation
        index=newpbin->add_ket(occ);
        newpbin->p[index]=newpbin->p[index]+p[i];
        delete occ;
    }

    newpbin->N=N;

    // Return new bin
    return newpbin;
}


//--------------------------------------------
//
// Adds Gaussian white noise to the stored statistics
//
//--------------------------------------------
p_bin *p_bin::white_noise(double stdev){
//  double stdev;       // Standard deviation of the noise
//  Variables
    double value;       // Sampled value of the normal distribution
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index


    newpbin=this->clone();
    for(i=0;i<newpbin->nket;i++){
        value=grand(0.0,stdev);
        newpbin->p[i]=newpbin->p[i]+newpbin->N*value;
        if (newpbin->p[i]<0) newpbin->p[i]=0.0;
    }

    // Return new bin
    return newpbin;
}


//--------------------------------------------
//
// Computes the effect of the detector blinking
//
//--------------------------------------------
p_bin *p_bin::blink(int S, qocircuit* qoc){
//  int S;              // Number of iterations to perform the computation
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    s;           // Wavepacket number
    int    blnk;        // Blink channel 0=No/1=yes
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    double prob;        // Blinking probability
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If the number of iterations is zero just return the input
    if(S==0){
        newpbin=this->clone();
        return newpbin;
    }

    // Calculate blinking
    newpbin=new p_bin(nph, nlevel,maxket,vis);
    for(j=0;j<S;j++){
        for(i=0;i<nket;i++){
            occ=new int[nlevel]();
            for(s=0;s<qoc->ns;s++){
            for(ch=0;ch<qoc->ndetc;ch++){
                // For each channel calculate the probability of blinking at some determined time
                prob=urand();
                if(prob>=qoc->det_par(0,ch)) blnk=0;
                else blnk=1;

                // Store the level occupation depending if the channel is blinking at that moment
                for(m=0;m<qoc->nm;m++){
                    l=qoc->i_idx[ch][m][s];
                    k=0;
                    while(vis[k]!=l) k=k+1;

                    if(blnk==0) occ[k]=occ[k]+ket[i][k];
                    else occ[k]=0;
                }
            }}

            //Loss channel do not blink. Their occupations are just copied
            for(s=0;s<qoc->ns;s++){
            for(ch=qoc->ndetc;ch<qoc->nch;ch++){
                for(m=0;m<qoc->nm;m++){
                    l=qoc->i_idx[ch][m][s];
                    k=0;
                    while(vis[k]!=l) k=k+1;

                    occ[k]=occ[k]+ket[i][k];
                }
            }}

            // Store level occupation
            index=newpbin->add_ket(occ);
            newpbin->p[index]=newpbin->p[index]+p[i];
            delete occ;
        }
    }
    newpbin->N=S*N;

    // Return new bin
    return newpbin;
}


//--------------------------------------------
//
// Computes the effect of the addition of dark counts
//
//--------------------------------------------
p_bin *p_bin::dark_counts(int S, qocircuit* qoc){
//  int S;              // Number of iterations to perform the computation
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    number;      // Number of photons created by dark counts
    int    total;       // Total number of photons create by dark counts
    int    nonzero;     // Number of channel with non-zero additions due to dark counts
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If the number of iterations is zero just return the input
    newpbin=this->clone();
    if(S==0) return newpbin;

    // Calculate Dark counts
    nonzero=0;
    for(i=0;i<newpbin->nket;i++) newpbin->p[i]=S*newpbin->p[i];
    for(i=0;i<S;i++){
        occ=new int[nlevel]();
        total=0;
        for(ch=0;ch<qoc->ndetc;ch++){
        for(m=0;m<qoc->nm;m++){
                // For each channels the number of extra dark counts photons is obtained
                number=prand(qoc->det_par(1,ch));
                l=qoc->i_idx[ch][m][0];
                k=0;
                while(vis[k]!=l) k=k+1;

                occ[k]=number;
                total=total+number;
        }}
        // If we add some new photon the new ket is stored
        if(total>0){
            index=newpbin->add_ket(occ);
            newpbin->p[index]=newpbin->p[index]+(double)N;
            nonzero=nonzero+1;
        }
        delete occ;
    }
    newpbin->N=newpbin->N*S;

    // Return new bin
    return newpbin;
}


//----------------------------------------
//
// Computes the effect of conditional detection
// from a table definition.
//
//----------------------------------------
p_bin *p_bin:: post_select_cond(int ndec, mati def,qocircuit *qoc){
//  int        ndec;                // Number of detector. ndec<nch because loss channels have no detector
//  mati       def;                 // Definition of the conditions of detection
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    int        first;               // Is this the first post-selection 1='Yes'/0='No'
    int        nph;                 // Number of photons
    int        eph;                 // Expected number of photons
    int        maxch;               // Larger channel number
    int        selbase;             // Base of the projector definition hash table
    int        nempty;              // Number of channels with zero photons
    int        nentry;              // Number of entries in the projector definition
    int        nprj;                // Number of projectors
    int        iph;                 // Index of photons
    int        ich;                 // Channel index
    int        im;                  // Mode index
    int        is;                  // "Time" index
    int       *tim;                 // "Time" configuration vector
    int       *pol;                 // Polarization configuration vector
    veci       ch;                  // List of channels by photons
    veci       pch;                 // List of polarization by photons
    hterm      select;              // Post-selection condition
    p_bin     *post_selected;       // New post-selected state for a single projector
    p_bin     *conditioned;         // New post-selected bin
    projector *prj;                 // Projector
//  Hash tableS
    int        key[3];              // Key of the projector definition
    int       *keyprj;              // Key of the projector
    long int   selvalue;            // Projector definition index
    long int   prjvalue;            // Projector index
    thash      selhash;             // Definitions hash table
    thash      prjhash;             // Projectors hash table.
    thash::const_iterator vselhash; // Row hash value.
    thash::const_iterator vprjhash; // Row hash value.
//  Auxiliary index
    int        i;                   // Aux index.
    int        j;                   // Aux index
    int        k;                   // Aux index
    int        l;                   // Aux index


    // Init variables and reserve memory.
    conditioned=new p_bin(this->nph, nlevel,maxket,this->vis);
    first=1;

    // Calculate:
    // How many photons.
    // How many empty channels
    // Maximum number of photons in channel
    nph=0;
    nempty=0;
    for(ich=0;ich<ndec;ich++){
        nph=nph+def(1,ich);
        maxch=max(maxch,def(0,ich));
        if(def(1,ich)==0) nempty=nempty+1;

    }
    selbase=max(maxch,max(qoc->nm,qoc->ns));

    // Create a list for every photon in which channel they are.
    // Empty channels that are post-selected to zero are appended at the end.
    ch.resize(nph+nempty);
    pch.resize(nph+nempty);
    k=0;
    l=0;
    for(ich=0;ich<ndec;ich++){
        for(iph=0;iph<def(1,ich);iph++){
            ch(k)=def(0,ich);
            pch(k)=def(2,ich);
            k++;
        }

        if(def(1,ich)==0){
            ch(nph+l)=def(0,ich);
            pch(nph+l)=-1;
            l++;
        }
    }

    // Initialize projector aux variables.
    nprj=0;
    keyprj=new int[ndec*qoc->nm*qoc->ns]();

    // From the detector definition of conditional detection create
    // all the projections necessary considering all the possible
    // polarizations and times.

    // LOOP Polarization¡
    pol=new int[nph+1]();
    while(pol[nph]==0){

        // LOOP "Time"
        tim=new int[nph+1]();
        while(tim[nph]==0){
            // Create Projector definition
            eph=0;
            nentry=0;
            selhash.clear();
            select.setZero(4,ndec*qoc->nm*qoc->ns);
            for(iph=0;iph<(nph+nempty);iph++){
            for(im=0;im<qoc->nm;im++){
            for(is=0;is<qoc->ns;is++){
                key[0]=ch(iph);
                key[1]=im;
                key[2]=is;


                // We have to check if for those photons the entry in select
                // already exists. If exist we update it otherwise it has to be
                // created.
                selvalue=hashval(key,3,selbase);
                vselhash=selhash.find(selvalue);
                if(vselhash==selhash.end()){
                    k=nentry;
                    selhash[selvalue]=nentry;
                    nentry=nentry+1;
                }else{
                    k=vselhash->second;
                }

                select(0,k)=ch(iph);
                select(1,k)=im;
                select(2,k)=is;
                if((iph<nph)&&(im==pol[iph])&&(is==tim[iph])&&((im==pch(iph))||(pch(iph)==-1))){
                    select(3,k)=select(3,k)+1;
                    eph=eph+1;
                }
            }}}


            // Create projector.
            // Note that polarization requirement may discard photons. In this case we do not
            // create the projector.
            if(eph==nph){
                // Check that we have not generated that projector before.
                // Photons are indistinguishable therefore in the way this loop
                // is built repetitions may appear. If it has not been generated
                // before the projector is created
                for(i=0;i<(ndec*qoc->nm*qoc->ns);i++){
                    keyprj[i]=select(3,i);
                }
                prjvalue=hashval(keyprj,ndec*qoc->nm*qoc->ns, this->nph);
                vprjhash=prjhash.find(prjvalue);
                // If the projector has not been created before
                if(vprjhash==prjhash.end()){
                    // We updated the entry in the corresponding hash table
                    prjhash[prjvalue]=nprj;
                    nprj=nprj+1;

                    // We create the projector
                    prj=new projector(this->nph, this->nlevel,2,vis);
                    prj->add_term(1.0,select,qoc);

                    // Apply the projector to the input state and store the result
                    post_selected=this->post_selection(prj);

                    if(first==1){
                        delete conditioned;
                        conditioned=new p_bin(post_selected->nph, post_selected->nlevel,maxket,post_selected->vis);
                        first=0;
                    }

                    conditioned->add_bin(post_selected);

                    // Free memory
                    delete post_selected;
                    delete prj;
                }
            }

            // ADVANCE "Time" in LOOP
            tim[0]=tim[0]+1;
            for(j=0;j<nph;j++){
                if(tim[j]>=qoc->ns){
                    tim[j]=0;
                    tim[j+1]=tim[j+1]+1;
                }
            }
        }

        // ADVANCE Polarization in LOOP
        pol[0]=pol[0]+1;
        for(j=0;j<nph;j++){
            if(pol[j]>=qoc->nm){
                pol[j]=0;
                pol[j+1]=pol[j+1]+1;
            }
        }
        // At the end of the loop free memory
        delete[] tim;
    }
    // At the end of the loop free memory
    delete[] pol;

    // Free memory;
    delete[] keyprj;
    selhash.clear();
    prjhash.clear();

    conditioned->N=N;

    // Return new bin
    return conditioned;
}


//----------------------------------------
//
// Computes the effect of conditional detection
// from the circuit stored definition.
//
//----------------------------------------
p_bin *p_bin:: compute_cond(qocircuit *qoc){
//  qocircuit *qoc;                 // Circuit where the detectors are placed
//  Variables
    p_bin     *conditioned;         // New post-selected bin


    conditioned=this->post_select_cond(qoc->ncond,qoc->det_def,qoc);
    return conditioned;
}


//----------------------------------------
//
// Traces out the channels given by a list.
//
//----------------------------------------
p_bin *p_bin:: remove_channels(veci ch, qocircuit *qoc){
//  veci       ch;                  // List of channels to be traced out
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    int        ich;                 // Index of the list ch
    int        nph;                 // Occupation of each level
    mati       select;              // Post-selection condition  for each level
    p_bin     *aux;                 // Auxiliary post_selected bin
    p_bin     *next;                // Auxiliary partially traced out bin
    p_bin     *removed;             // Final traced out bin


    // Init variables
    select.resize(3,1);
    removed=this->clone();
    next=new p_bin(removed->nph, removed->nlevel,maxket,removed->vis); // To avoid warning. Not really needed

    // For each channel
    for(ich=0;ich<ch.size();ich++){
        // Post-select each of the photon occupations. Trace out
        for(nph=0;nph<=this->nph;nph++){
            select(0,0)=ch(ich);
            select(1,0)=nph;
            select(2,0)=-1;
            aux=removed->post_select_cond(select.cols(),select,qoc);

            if(nph==0) {
                delete next;
                next=new p_bin(aux->nph, aux->nlevel,maxket,aux->vis);
            }
            next->add_bin(aux);
            delete aux;
        }
        // Update total bin
        delete removed;
        next->N=N;
        removed=next->clone();
    }

    // Free memory
    delete next;

    // Return new bin
    return removed;

}


//----------------------------------------
//
// Computes the effect of detection considering losses
//
//----------------------------------------
p_bin *p_bin:: compute_loss(qocircuit *qoc){
//  qocircuit *qoc;              // Circuit where the detector is placed
//  Variables
    int        fchloss;          // First loss channel
    int        nchloss;          // Number of loss channels
    veci       chloss;           // List of loss channels
    p_bin     *lossed;           // New post-selected bin
//  Auxiliary index
    int        i;                // Aux index


    // Init variables
    nchloss=qoc->nch/2;
    fchloss=qoc->nch/2;

    // Create list of loss channels
    chloss.resize(nchloss);
    for(i=0;i<nchloss;i++) chloss(i)=fchloss+i;

    // Remove the loss channels.
    lossed=this->remove_channels(chloss,qoc);

    // Return resulting pbin
    return lossed;
}


//----------------------------------------
//
// Computes the effect of ignoring channels
//
//----------------------------------------
p_bin *p_bin:: compute_ignored(qocircuit *qoc){
//  qocircuit *qoc;              // Circuit where the detector is placed
//  Variables
    int    nignored;             // Number of ignored channels
    veci   chign;                // List of ignored channels
    p_bin *ignored;              // New post-selected bin
//  Auxiliary index
    int i;                       // Aux index


    // Init variables
    nignored=qoc->nignored;

    // Create list of ignored channels
    chign.resize(nignored);
    for(i=0;i<nignored;i++) chign(i)=qoc->ch_ignored(i);

    // Remove the ignored channels.
    ignored=this->remove_channels(chign,qoc);

    // Return resulting pbin
    return ignored;
}


//----------------------------------------
//
// Remove frequencies and leaves only times
//
//----------------------------------------
p_bin *p_bin::remove_freq(qocircuit* qoc){
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    is;          // Packet number
    int    it;          // Packet number
    int    ip;          // Period number
    int    nt;          // Number of times
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If for some reason there are not packets defined
    // then the input is returned
    if(qoc->emitted->pack_def.cols()==0){
        newpbin=this->clone();
        return newpbin;
    }

    // Set up number of times.
    nt=qoc->emitted->times.size();

    // Calculate counts in a channel-polarization
    // VISIBILITY information is to remember the post-selected channels.
    newpbin=new p_bin(nph, nlevel,nket,vis);
    for(i=0;i<nket;i++){
        occ=new int[qoc->nlevel]();
        for(j=0;j<nlevel;j++){
            ch=qoc->idx[vis[j]].ch;
            m=qoc->idx[vis[j]].m;
            is=qoc->idx[vis[j]].s%qoc->nsp;
            ip=qoc->idx[vis[j]].s/qoc->nsp;

            if(ket[i][j]>0){ // We prevent to check levels that we don't have packet defined when we reserve more packets than we use
                // Redefine the index as time instead of packets
                it=qoc->emitted->pack_def(0,is)+ip*nt;
                l=qoc->i_idx[ch][m][it];

                // Store
                k=0;
                while(vis[k]!=l) k=k+1;
                occ[k]=occ[k]+ket[i][j];
            }
        }
        index=newpbin->add_ket(occ);
        newpbin->p[index]=newpbin->p[index]+p[i];
        delete occ;
    }
    newpbin->N=N;

    // Return new bin.
    return newpbin;
}


//----------------------------------------
//
// Classify period
//
//----------------------------------------
p_bin *p_bin::classify_period(qocircuit* qoc){
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    ip;          // Packet period
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // If for some reason there are not packets defined
    // then the input is returned
    if(qoc->emitted->pack_def.cols()==0){
        newpbin=this->clone();
        return newpbin;
    }

    // Calculate counts in a channel-polarization
    // VISIBILITY information is to remember the post-selected channels.
    newpbin=new p_bin(nph,nlevel,nket,vis);
    for(i=0;i<nket;i++){
        occ=new int[qoc->nlevel]();
        for(j=0;j<nlevel;j++){
            ch=qoc->idx[vis[j]].ch;
            m=qoc->idx[vis[j]].m;
            ip=qoc->idx[vis[j]].s/qoc->nsp;

            if(ket[i][j]>0){ // We prevent to check levels that we don't have packet defined when we reserve more packets than we use
                // Redefine the index as period instead of packets
                l=qoc->i_idx[ch][m][ip];

                // Store
                k=0;
                while(vis[k]!=l) k=k+1;
                occ[k]=occ[k]+ket[i][j];
            }
        }
        index=newpbin->add_ket(occ);
        newpbin->p[index]=newpbin->p[index]+p[i];
        delete occ;
    }
    newpbin->N=N;

    // Return new bin.
    return newpbin;
}


//----------------------------------------
//
// Performs the count of photons in a channel
// independently of time and frequency
//
//----------------------------------------
p_bin *p_bin::perform_count(qocircuit* qoc){
//  qocircuit *qoc;     // Circuit where the detectors are defined
//  Variables
    int    ch;          // Channel number
    int    m;           // Mode number
    int    index;       // Index of a ket in this set of bins
    int   *occ;         // Level occupation
    p_bin *newpbin;     // New output probability bin
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int    l;           // Aux index


    // Calculate counts in a channel-polarization
    newpbin=new p_bin(nph,nlevel,nket,vis);
    for(i=0;i<nket;i++){
        occ=new int[nlevel]();
        for(j=0;j<nlevel;j++){
            // We put all the counts in
            // the same time/wavepacket
            ch=qoc->idx[vis[j]].ch;
            m=qoc->idx[vis[j]].m;
            l=qoc->i_idx[ch][m][0];

            // Store
            k=0;
            while(vis[k]!=l) k=k+1;
            occ[k]=occ[k]+ket[i][j];
        }
        index=newpbin->add_ket(occ);
        newpbin->p[index]=newpbin->p[index]+p[i];
        delete occ;
    }
    newpbin->N=N;

    // Return new bin.
    return newpbin;
}



//----------------------------------------
//
// Remove all the packet degrees of freedom
// except the 0.
//
//----------------------------------------
p_bin *p_bin:: remove_time(qocircuit *qoc){
//  qocircuit *qoc;              // Circuit where the detector is placed
//  Variables
    int    index;                // Index of a ket in this set of bins
    int    newnlevel;            // New number of levels
    int   *auxket;               // Auxiliary ket
    int   *newvis;               // New visibility vector
    int   *isincluded;           // Levels included
    p_bin *auxlist;              // Auxiliary ket list (to be returned)
//  Auxiliary index
    int    i;                    // Aux index
    int    j;                    // Aux index
    int    k;                    // Aux index


    // Init variables
    newnlevel=0;
    newvis=new int[nlevel]();
    isincluded=new int[nlevel]();
    for(i=0;i<nlevel;i++){
        if(qoc->idx[vis[i]].s==0){
                isincluded[i]=1;
                newvis[newnlevel]=vis[i];
                newnlevel=newnlevel+1;
        }
    }
    auxlist=new p_bin(nph,newnlevel,maxket,newvis);
    delete[] newvis;

    // Perform operations
    auxket=new int[newnlevel]();
    for(i=0;i<nket;i++){
        k=0;
        for(j=0;j<nlevel;j++){
            if(isincluded[j]==1){
                auxket[k]=ket[i][j];
                k++;
            }
        }
        index=auxlist->add_ket(auxket);
        auxlist->p[index]=auxlist->p[index]+p[i];
    }
    auxlist->N=N;

    // Free memory
    delete auxket;
    delete[] isincluded;

    // Return new bin
    return auxlist;
}


//----------------------------------------
//
// Probability value a bin from it definition
// It uses qodev instead of qocircuit
//
//----------------------------------------
double p_bin::prob(mati def,qodev *dev){
//  mati       def;          // Definition of a probability bin in human readable form
//  qodev     *dev;          // Device where the detector is placed


    return this->prob(def, dev->circ);
}


//---------------------------------------------------------------------------
//
//  Encode a p_bin into a qubit representation using path encoding
//  Circuit version.
//
//----------------------------------------------------------------------------
p_bin *p_bin::translate(mati qdef,qocircuit *qoc){
// mati       qdef; // Integer matrix with the qubit to channel definitions
// qocircuit *qoc;  // Quatum optical circuit to which this p_bin is referred
// Variables
    bool   valid;   // Is the ket valid for codification true=Yes/false=No
    bool   print;   // A warning has been printed. True=Yes/False=No
    int    val0;    // value of channel 0
    int    val1;    // value of channel 1
    int    qval;    // qubit equivalent value of channels 0 and 1.
    int    pos;     // Position of the new stored value
    int    nvalid;  // Number of encoded bins
    int   *values;  // Vector of the values of each qubit
    p_bin *qbin;    // qubit encoded probability bin
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index
    int    l;       // Aux index
    int    m;       // Aux index
    int    n;       // Aux index


    // Reserve and initialize memory
    qbin=new p_bin(1,qdef.cols(),maxket);

    // Check translation conditions
    if((qoc->nm>1)||(qoc->ns>1)){
        cout << "Translate error: The circuit number of modes nm and packets ns can only be one for encoding." << endl;
        return qbin;
    }

    // Encode each bin
    nvalid=0;
    print=false;
    for(i=0;i<nket;i++){
        values=new int[qdef.cols()]();
        valid=true;
        for(j=0;j<qdef.cols();j++){
            // Find the channels
            m=qoc->i_idx[qdef(0,j)][0][0];
            n=qoc->i_idx[qdef(1,j)][0][0];

            // Read the values
            // qdef is given over circuit definition.
            // Levels number may change after post-selection
            k=0;
            l=0;
            while((vis[k]!=m)&&(k<nlevel)) k=k+1; // Not efficient but small search
            while((vis[l]!=n)&&(l<nlevel)) l=l+1; // Not efficient but small search
            val0=ket[i][k];
            val1=ket[i][l];

            // Encode the values
            qval=-1;
            if((val0==0)&&(val1==1)) qval=0;
            if((val0==1)&&(val1==0)) qval=1;
            if(qval<0) valid=false;
            values[j]=qval;
        }

        if(valid==true){
            pos=qbin->add_count(values);
            qbin->p[pos]=p[i];
            if((pos!=nvalid)&&(print==false)){
                cout << "Translate: Warning! encoding leads to collision. Invalid result" << endl;
                print=true;
            }
            nvalid=nvalid+1;
        }
        delete values;
    }

    // Normalize bins
    qbin->N=N;

    // Return encoded bins
    return qbin;
}


//---------------------------------------------------------------------------
//
//  Encode a p_bin into a qubit representation using path encoding.
//  Device version.
//
//----------------------------------------------------------------------------
p_bin *p_bin::translate(mati qdef,qodev *dev){
//  mati       qdef; // Integer matrix with the qubit to channel definitions
//  qodev     *dev;  // Quatum optical device to which this p_bin is referred


    return translate(qdef,dev->circ);
}


//---------------------------------------------------------------------------
//
//  Encode a p_bin into a qubit representation using polarization encoding.
//  Circuit version.
//
//----------------------------------------------------------------------------
p_bin *p_bin::pol_translate(veci qdef,qocircuit *qoc){
// veci       qdef; // Integer vector with the qubit to channel definitions
// qocircuit *qoc;  // Quatum optical circuit to which this p_bin is referred
// Variables
    bool   valid;   // Is the ket valid for codification true=Yes/false=No
    bool   print;   // A warning has been printed. True=Yes/False=No
    int    valH;    // value of channel 0
    int    valV;    // value of channel 1
    int    qval;    // qubit equivalent value of channels 0 and 1.
    int    pos;     // Position of the new stored value
    int    nvalid;  // Number of encoded bins
    int   *values;  // Vector of the values of each qubit
    p_bin *qbin;    // qubit encoded probability bin
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index
    int    l;       // Aux index
    int    m;       // Aux index
    int    n;       // Aux index


    // Reserve and initialize memory
    qbin=new p_bin(1,qdef.size(),maxket);

    // Check translation conditions
    if((qoc->nm!=2)||(qoc->ns>1)){
        cout << "Pol_translate error: The circuit number of modes nm and packets ns should be nm=2 and ns=1 for encoding." << endl;
        return qbin;
    }

    // Encode each bin
    nvalid=0;
    print=false;
    for(i=0;i<nket;i++){
        values=new int[qdef.size()]();
        valid=true;
        for(j=0;j<qdef.size();j++){
            // Find the channels
            m=qoc->i_idx[qdef(j)][H][0];
            n=qoc->i_idx[qdef(j)][V][0];

            // Read the values
            // qdef is given over circuit definition.
            // Levels number may change after post-selection
            k=0;
            l=0;
            while((vis[k]!=m)&&(k<nlevel)) k=k+1; // Not efficient but small search
            while((vis[l]!=n)&&(l<nlevel)) l=l+1; // Not efficient but small search
            valH=ket[i][k];
            valV=ket[i][l];

            // Encode the values
            qval=-1;
            if((valH==1)&&(valV==0)) qval=0;
            if((valH==0)&&(valV==1)) qval=1;
            if(qval<0) valid=false;
            values[j]=qval;
        }

        if(valid==true){
            pos=qbin->add_count(values);
            qbin->p[pos]=p[i];
            if((pos!=nvalid)&&(print==false)){
                cout << "Pol_translate: Warning! encoding leads to collision. Invalid result" << endl;
                print=true;
            }
            nvalid=nvalid+1;
        }
        delete values;
    }

    // Normalize bins
    qbin->N=N;

    // Return encoded bins
    return qbin;
}


//---------------------------------------------------------------------------
//
//  Encode a p_bin into a qubit representation using polarization encoding.
//  Device version.
//
//----------------------------------------------------------------------------
p_bin *p_bin::pol_translate(veci qdef,qodev *dev){
//  veci       qdef; // Integer vector with the qubit to channel definitions
//  qodev     *dev;  // Quatum optical device to which this p_bin is referred


    return pol_translate(qdef,dev->circ);
}
