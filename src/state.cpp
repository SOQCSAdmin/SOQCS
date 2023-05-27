//======================================================================================================
// File state.cpp
//
// STATE DEFINITION LIBRARY.
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "state.h"


//----------------------------------------
//
//  Creates a ket list.
//  The maximum number of kets is set by default.
//
//----------------------------------------
ket_list::ket_list(int i_level){
//  int i_level     // Number of levels to describe a ket


    create_ket_list(i_level,DEFSTATEDIM);
}


//-----------------------------------------------------
//
//  Creates a ket list specifying the maximum number of kets.
//
//-----------------------------------------------------
ket_list::ket_list(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe a ket
//  int i_maxket    // Maximum number of different kets in the list


    create_ket_list(i_level,i_maxket);
}


//-----------------------------------------------------
//
//  Creates a ket list specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
ket_list::ket_list(int i_level, int i_maxket, int *i_vis){
//  int  i_level     // Number of levels to describe a ket
//  int  i_maxket    // Maximum number of different kets in the list
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels
//  Auxiliary  index
    int  g;           // Aux index


    create_ket_list(i_level,i_maxket);
    for(g=0;g<i_level;g++) vis[g]=i_vis[g];
}


//----------------------------------------
//
//  Auxiliary private method to create a ket list
//  Each ket is distinguished by level occupations
//
//----------------------------------------
void ket_list::create_ket_list(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe a ket
//  int i_maxket    // Maximum number of different kets in the list
//  Auxiliary index
    int i;          // Aux index


    // Init ket information
    nket=0;
    nlevel=i_level;
    maxket=i_maxket;

    // Create and int amplitude/coefficient and ket structures
    ket=new int*[maxket];
    for(i=0;i<maxket;i++){
        ket[i]=new int[nlevel]();
    }

    // Compute the trivial print visibility vector.
    vis=new int[nlevel];
    for(i=0;i<nlevel;i++){
        vis[i]=i;
    }
}


//----------------------------------------
//
//  Destroys a ket_list
//
//----------------------------------------
ket_list::~ket_list(){
//  Auxiliary index
    int i;          // Aux index


    // Liberate memory of amplitude/coefficients and ket structures
    for(i=0;i<maxket;i++){
        delete[] ket[i];
    }

    //Free memory
    delete[] ket;
    delete[] vis;
    // Clear hash table
    ketindex.clear();
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
ket_list *ket_list::clone(){
    // Variable
    ket_list *aux;  // Auxiliary state recipient of the copy
    // Auxiliary index
    int i;          // Aux index
    int j;          // Aux index


    aux=new ket_list(nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }
    aux->ketindex=ketindex;
    return aux;
}


//----------------------------------------
//
// Finds the position of a ket in the list
//
//----------------------------------------
int ket_list::find_ket(int *occ){
//  int *occ;                    // Occupation of those levels
//  Variables
    long long int value;         // Hash value in a dynamic base dictionary
    thash::const_iterator hashv; // Hash value constant iterator.


    // Update amplitude and occupation
    // Store
    value=hashval(occ,nlevel);
    hashv=ketindex.find(value);
    // If the ket does not exist yet on the output create new one.
    if(hashv==ketindex.end()){
        return -1;
    }else{
        return hashv->second;
    }
}


//----------------------------------------
//
// Finds the position of a ket in the list
// The ket is defined in human readable form
//
//----------------------------------------
int ket_list::find_ket(mati def,qocircuit *qoc){
//  mati       def;            // Definition of the ket that is searched
//  qocircuit *qoc;            // Circuit to which this list is referred
//  Variables
    ket_list  *aux;            // Bra state created from the definition def
//  Dictionary variables
    int        index;          // Index of rows


    // Create bra reference ket
    aux= new ket_list(nlevel,1,vis);
    aux->add_ket(def,qoc);

    // Search the reference
    index=find_ket(aux->ket[0]);

    // Free memory
    delete aux;

    // Return value
    return index;
}


//----------------------------------------
//
//  Adds a new ket to the list state
//
//----------------------------------------
int ket_list::add_ket(int *occ){
//  int *occ;                    // Level occupation
//  Variables
    int  index;                  // Index/List position to be returned
    long long int value;         // Hash value in a dynamic base dictionary
    thash::const_iterator hashv; // Hash value constant iterator.
//  Auxiliary index
    int  i;                      // Aux index


    // Check and warn about memory limits.
    if(nket>=maxket){
        cout << "Ket list add ket error #1: Memory limit exceeded." << endl;
        return -1;
    }

    // Update amplitude and occupation
    //Store
    value=hashval(occ,nlevel);
    hashv=ketindex.find(value);

    if(hashv==ketindex.end()){
        index=nket;
        ketindex[value]=nket;
        for(i=0;i<nlevel;i++){
            ket[nket][i]=occ[i];

        }

        nket=nket+1;
    }else{
        index=hashv->second;
    }

    // Return storage index
    return index;
}


//---------------------------------------------------
//
// Adds a new ket to the list (human readable form)
//
//---------------------------------------------------
int ket_list::add_ket(hterm term, qocircuit *qoc){
//  hterm term                   // Matrix with the state specification in human readable form
//  qocircuit *qoc               // Circuit to which the the term is referred (to translate from human to numeric form)
//  Variables
    int  index;                  // Index to be returned
    int  olevel;                 // Level index
    int  ilevel;                 // Level index
    int *occ;                    // Occupation
    int *ivis;                   // Inverse visibility.
//  Auxiliary index
    int  i;                      // Aux index
    int  num;                    // Number of columns in term


    // Check and warn about memory limits.
    if(nket>=maxket){
        cout << "Ket list add ket error #2: Memory limit exceeded." << endl;
        exit(0);
    }

    // Initializations
    ivis=new int[qoc->nlevel];
    for(i=0;i<qoc->nlevel;i++) ivis[i]=-1;
    for(i=0;i<nlevel;i++) ivis[vis[i]]=i;
    occ=new int[nlevel]();
    for(i=0;i<nlevel;i++) occ[i]=ket[maxket-1][0];

    num=term.cols();
    // Update occupation
    for(i=0;i<num;i++){
        // The number of rows of the matrix determines its format.
        switch (term.rows()){
            case 4:     //We fully describe the state.
                if (term(0,i)>=qoc->nch){ cout  << "add_ket: Warning! channel not defined." << endl; return -1;}
                if (term(1,i)>=qoc->nm ){ cout  << "add_ket: Warning! mode/polarization not defined." << endl; return -1;}
                if (term(2,i)>=qoc->ns ){ cout  << "add_ket: Warning! packet not defined."  << endl; return -1;}

                olevel=qoc->i_idx[term(0,i)][term(1,i)][term(2,i)];
                ilevel=ivis[olevel];
                if(ilevel<0){
                    cout << "add_ket error: That level no longer exist in this Hilbert space" << endl;
                    exit(0);
                }
                occ[ilevel]=term(3,i);
                break;
            case 3:     // We skip bunch (time and frequency)
                        // There are relevant in just a few problems.
                if (term(0,i)>=qoc->nch){ cout  << "add_ket: Warning! channel not defined." << endl; return -1;}
                if (term(1,i)>=qoc->nm ){ cout  << "add_ket: Warning! mode/polarization not defined." << endl; return -1;}

                olevel=qoc->i_idx[term(0,i)][term(1,i)][0];
                ilevel=ivis[olevel];
                if(ilevel<0){
                    cout << "add_ket error: That level no longer exist in this Hilbert space" << endl;
                    exit(0);
                }
                occ[ilevel]=term(2,i);
                break;
            case 2:     // We consider there is no polarization in this problem. (Just channel and occupation).
                if (term(0,i)>=qoc->nch){ cout  << "add_ket: Warning! channel not defined." << endl; return -1;}

                olevel=qoc->i_idx[term(0,i)][0][0];
                ilevel=ivis[olevel];
                if(ilevel<0){
                    cout << "add_ket error: That level no longer exist in this Hilbert space" << endl;
                    exit(0);
                }
                occ[ilevel]=term(1,i);
                break;
            case 1:     //  Like before we can just input occupations.
                        //  However in this case we must define all of them.
                if (num>qoc->nlevel){ cout << "Add ket: Warning! level not defined." << endl; return -1;}
                olevel=qoc->i_idx[i][0][0];
                ilevel=ivis[olevel];
                if(ilevel<0){
                    cout << "add_ket error: That level no longer exist in this Hilbert space" << endl;
                    exit(0);
                }
                occ[ilevel]=term(0,i);
                break;
            default:    // Something has gone wrong.
                break;
        }
    }

    // Update list
    index=add_ket(occ);

    // Free memory
    delete []occ;
    delete []ivis;

    // Return index.
    return index;
}


//-------------------------------------------------------
//
// Remove levels with wave-packet indexes greater than zero.
// Needed for some specific operations
//
//-------------------------------------------------------
ket_list *ket_list:: remove_time(qocircuit *qoc){
//  qocircuit *qoc;                 // Circuit where that is referred by this list
//  Variables
    int  newnlevel;                 // New number of levels
    int *auxket;                    // Auxiliary ket
    int *newvis;                    // New visibility vector
    int *isincluded;                // Levels included
    ket_list *auxlist;              // Auxiliary ket list (to be returned)
//  Auxiliary index
    int i;                          // Aux index
    int j;                          // Aux index
    int k;                          // Aux index


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

    auxlist=new ket_list(newnlevel,maxket,newvis);
    delete[] newvis;

    // Compute operation
    auxket=new int[newnlevel]();
    for(i=0;i<nket;i++){
        k=0;
        for(j=0;j<nlevel;j++){
            if(isincluded[j]==1){
                auxket[k]=ket[i][j];
                k++;
            }
        }
        auxlist->add_ket(auxket);
    }

    // Free memory
    delete auxket;
    delete[] isincluded;

    // Return value
    return auxlist;
}


//-------------------------------------------------------
//
//  Clear the list (without destroying it).
//
//-------------------------------------------------------
void ket_list::clear_kets(){


    nket=0;
    ketindex.clear();
}


//--------------------------------------------
//
//  Prints a ket in a default format
//
//--------------------------------------------
void  ket_list::prnt_ket(int iket){
//  int iket;           // Ket number to be printed


    prnt_ket(iket, DEFFORMAT, false, nullptr);
}


//--------------------------------------------
//
//  Prints a ket
//
//--------------------------------------------
void  ket_list::prnt_ket(int iket, int format, qocircuit *qoc){
//  int iket;           // Ket number to be printed
//  int format;         // Integer that defines the notation used when printing the ket
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)


    prnt_ket(iket, format, false, qoc);
}


//--------------------------------------------
//
//  Prints a ket
//
//--------------------------------------------
void  ket_list::prnt_ket(int iket, int format, bool loss, qocircuit *qoc){
//  int iket;           // Ket number to be printed
//  int format;         // Integer that defines the notation used when printing the ket
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  Variables
    int writeprev;      // Have we print a level of this ket previously? 1=Yes/0=No
    int lev;            // Level number as defined in circuit base. We do not rely in order anymore because post-selection may change it.
    int nchm;           // Half the number of channels.
//  Auxiliary index
    int k;              // Aux index


    cout  << CYAN << " | ";
    writeprev=0;
    for(k=0;k<nlevel;k++){
        lev=vis[k];
        if(qoc > (qocircuit *)nullptr){
            nchm=qoc->nch/2;
            // We print in different formats depending on the circuit flag for printing
            switch (format){
                case 0: // Straightforward form
                    if(writeprev==1) cout << ", ";
                    if((qoc->idx[lev].ch>=nchm)&&(loss)) cout  << BLUE;
                    if(ket[iket][k]>=0){
                        cout << ket[iket][k];
                    }else{
                        cout << "X";
                    }
                    cout  << CYAN;
                    writeprev=1;
                    break;


                case 1: // Human readable condensed form. A condensed form particularly apt for Baso Baset problem and similars.
                    if(ket[iket][k]>0){
                        if(writeprev==1) cout << ", ";
                        if((qoc->idx[lev].ch>=nchm)&&(loss)) cout  << RED;
                        if(ket[iket][k]>1) cout << "[" << ket[iket][k] << "]";
                        if(qoc->nm>1)    cout << pl[qoc->idx[lev].m];
                        if(qoc->ns>1)    cout << "(" << qoc->idx[lev].s << ")";
                        cout << qoc->idx[lev].ch;
                        cout  << CYAN;
                        writeprev=1;
                    }
                    break;
                default:
                    cout << "Prnt_state error: Format "<< format << " does not exist." << endl;
                    exit(0);
                    break;
            }
        }else{
            if(writeprev==1) cout << ", ";
            if(ket[iket][k]>=0){
                cout << ket[iket][k];
            }else{
                cout << "X";
            }
            writeprev=1;
        }
    }
    cout << " >"<< RESET;
}


//----------------------------------------
//
//  Creates a state
//  The maximum number of kets is set by default.
//
//----------------------------------------
state::state(int i_level):ket_list(i_level){
//  int i_level;     // Number of levels to describe the state


    ampl=new cmplx[maxket]();
}


//-----------------------------------------------------
//
//  Creates a state specifying the maximum number of kets.
//
//-----------------------------------------------------
state::state(int i_level, int i_maxket):ket_list(i_level,i_maxket){
//  int i_level     // Number of levels to describe the state
//  int i_maxket    // Maximum number of different kets in the summation


    ampl=new cmplx[maxket]();
}


//-----------------------------------------------------
//
//  Creates a state specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
state::state(int i_level, int i_maxket, int *i_vis):ket_list(i_level,i_maxket, i_vis){
//  int  i_level     // Number of levels to describe the state
//  int  i_maxket    // Maximum number of different kets in the state
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels


    ampl=new cmplx[maxket]();
}


//----------------------------------------
//
//  Destroys a state
//
//----------------------------------------
state::~state(){


    // Liberate memory of amplitude
    delete[] ampl;
}


//----------------------------------------
//
//  Copy a state
//
//----------------------------------------
state *state::clone(){
    // Variable
    state *aux;  // Auxiliary state
    // Auxiliary index
    int    i;    // Aux index
    int    j;    // Aux index


    aux=new state(nlevel,maxket);
    aux->nket=nket;
    for(i=0;i<nket;i++){
        aux->ampl[i]=ampl[i];
        for(j=0;j<nlevel;j++){
            aux->ket[i][j]=ket[i][j];
        }
    }
    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }

    aux->ketindex=ketindex;
    return aux;
}


//----------------------------------------
//
//  Clear a state
//
//----------------------------------------
void state::clear(){


    delete[] ampl;
    ampl=new cmplx[maxket]();
    clear_kets();
}


//----------------------------------------
//
//  Adds a new term to a state
//
//----------------------------------------
int state::add_term(cmplx i_ampl, int *occ){
//  cmplx i_ampl;       // Amplitude of the new term
//  int  *occ;          // Occupation of those levels
//  Variables
    int   index;        // Index to be returned


    index=add_ket(occ);
    if(index>=0) ampl[index]=ampl[index]+i_ampl;

    return index;
}


//---------------------------------------------------
//
// Adds a new term to a state ( human readable form)
//
//---------------------------------------------------
int state::add_term(cmplx i_ampl,hterm term, qocircuit *qoc){
//  cmplx i_ampl;                // Amplitude of the new term
//  hterm term                   // Matrix with the state specification in human readable form
//  qocircuit *qoc               // Circuit to which the the term is referred (to translate from human to numeric form)
//  Variables
    int   index;                 // Index to be returned


    index=add_ket(term,qoc);
    if(index>=0) ampl[index]=ampl[index]+i_ampl;

    return index;
}


//---------------------------------------------------
//
// Calculates the braket <prj|state>
//
//---------------------------------------------------
cmplx state::braket(state *bra){
//  state *bra;     // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    cmplx  rampl;   // Amplitude to be returned
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index


    rampl=0.0;
    // For each projector term
    for(i=0;i<bra->nket;i++){
            // Find the corresponding ket in the state
            j=find_ket(bra->ket[i]);
            if(j>=0){
                rampl=rampl+conj(bra->ampl[i])*ampl[j];
            }
    }

    // Return state.
    return rampl;
}


//---------------------------------------------------
//
// Normalizes the state
//
//---------------------------------------------------
void state::normalize(){
//  Variables
    double tot;     // Total/ Norm of the state
//  Auxiliary index
    int    i;       // Aux index


    tot=0.0;
    for(i=0;i<nket;i++){
        tot=tot+real(conj(ampl[i])*ampl[i]);
    }
    if(abs(tot)>xcut){
        for(i=0;i<nket;i++){
            ampl[i]=ampl[i]/sqrt(tot);
        }
    }

}


//---------------------------------------------------
//
// Changes the global phase of the state.
// Def defines the term we want to rephase to a real value.
//
//---------------------------------------------------
int state::rephase(hterm def,qocircuit *qoc){
//  hterm def;      // State definition
//  qocircuit *qoc; // Circuit to which def is referred
//  Variables
    double re;      // Real part
    double im;      // Imaginary part
    double phase;   // Phase.
    cmplx ampl_def; // Amplitude of def state
    int idx;        // Index. Position of the state on the ket list.
//  Auxiliary index
    int    i;       // Aux index


    idx=find_ket(def,qoc);
    if(idx<0) return idx;

    ampl_def=ampl[idx];
    re=real(ampl_def);
    im=imag(ampl_def);
    phase=atan2(im,re);

    for(i=0;i<nket;i++) ampl[i]=ampl[i]*exp(-jm*phase);
    return idx;
}


//--------------------------------------------
//
//  Post selection process of the output state
//
//--------------------------------------------
state *state::post_selection(state *prj){
//  state *prj;      // Projector with the description of the levels (and their occupations) to be post-selected
//  Variables
    int    npost;    // Number of levels in which post selection is performed
    int    selected; // Is this state selected 0=No/1=Yes
    int   *islincl;  // Is level included? 0=No/1=Yes
    int   *occ;      // Occupation
    state *nstate;   // New post-selected state
// Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index


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
    nstate=new state(nlevel-npost,maxket);
    k=0;
    for(l=0;l<nlevel;l++){
        if(islincl[l]==1){
            nstate->vis[k]=vis[l];
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
                nstate->add_term(ampl[j]*conj(prj->ampl[i]),occ);
            }
        }
    }

    // Free memory
    delete[] occ;
    delete[] islincl;

    // Return state.
    return nstate;
}


//--------------------------------------------
//
//  Post selection process of the output state
//  It ignores photons out of the detection window.
//  Intended for internal use only by dmatrix.
//
//--------------------------------------------
state *state::post_selection(state *prj,qocircuit *qoc){
//  state *prj;      // Projector with the description of the levels (and their occupations) to be post-selected
//  qocircuit *qoc;  // Circuit to which the projector is referred
//  Variables
    int    npost;    // Number of levels in which post selection is performed
    int    selected; // Is this state selected 0=No/1=Yes
    int    nwi;      // Initial period of the window of detection
    int    nwf;      // Last period (+1) of the window of detection
    int    ch;       // Channel
    int    s;        // Packet number
    int   *islincl;  // Is level included? 0=No/1=Yes
    int   *occ;      // Occupation
    state *nstate;   // New post-selected state
// Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index


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
    nstate=new state(nlevel-npost,maxket);
    k=0;
    for(l=0;l<nlevel;l++){
        if(islincl[l]==1){
            nstate->vis[k]=vis[l];
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
                ch=qoc->idx[k].ch;
                s=qoc->idx[k].s;
                if((qoc->losses==0)||(ch<qoc->nch/2)){
                    nwi=qoc->det_win(0,ch);
                    nwf=qoc->det_win(1,ch)+1;
                    if(qoc->det_win(0,ch)<0) nwi=0;
                    if(qoc->det_win(1,ch)<0) nwf=qoc->np+1;
                }else{
                    nwi=0;
                    nwf=qoc->np+1;
                }

                if ((ket[j][k]!=prj->ket[i][k])&&(prj->ket[i][k]>=0)&&(s>=nwi*qoc->nsp)&&(s<nwf*qoc->nsp))  selected=0;
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
                nstate->add_term(ampl[j]*conj(prj->ampl[i]),occ);
            }
        }
    }

    // Free memory
    delete[] occ;
    delete[] islincl;

    // Return state.
    return nstate;
}


//-------------------------------------------------------
//
// Post select a state using a "projector" initialized to zero
// in the corresponding channels defined by ch
//
//-------------------------------------------------------
state *state:: remove_empty_channels(veci ch, int it, qocircuit *qoc){
//  veci       ch;                // Channel list where the post-selection is performed
//  int        it;                // Allowed "times" 0=Only 0/ 1=All
//  qocircuit *qoc;               // Circuit where the detector is placed
//  Variables
    int        nl;                // List length
    int        nch;               // Detector condition definition number of columns to read
    int        ich;               // Channel index
    int        im;                // Mode index
    int        is;                // Packet index
    int       *ket;               // Ket
    hterm      select;            // Post-selection condition
    veci       llist;             // Level list
    projector *auxprj;            // Projector
    state     *auxstate;          // Auxiliary state
    state     *newstate;          // New post-selected state
//  Auxiliary index
    int        i;                 // Aux index
    int        l;                 // Aux index
    int        k;                 // Aux index


    // Init variables
    nch=ch.size();

    // Post select-channels
    k=0;
    if(nch>0) {
        select.setZero(4,nch*qoc->nm*qoc->ns);
        for(ich=0;ich<nch;ich++){
        for(im=0;im<qoc->nm;im++){
        for(is=0;is<qoc->ns;is++){

            select(0,k)=ch(ich);
            select(1,k)=im;
            select(2,k)=is;
            select(3,k)=0;
            k++;

        }}}

        auxprj=new projector(qoc->num_levels(),1,vis);
        auxprj->add_term(1.0,select,qoc);
        auxstate=this->post_selection(auxprj);
        delete auxprj;
    }else{
        auxstate=this->clone();
    }


    // Post-select "Times"
    if(it==0){
        nl=0;
        llist.resize(auxstate->nlevel);
        for(i=0;i<auxstate->nlevel;i++){
            l=auxstate->vis[i];
            if((qoc->idx[l].s)>0){
                llist(nl)=i;
                nl=nl+1;
            }
        }

        ket= new int[auxstate->nlevel];

        for(i=0;i<auxstate->nlevel;i++) ket[i]=-1;
        for(i=0;i<nl;i++)  ket[llist(i)]=0;

        auxprj=new projector(auxstate->nlevel,1,auxstate->vis);
        auxprj->add_term(1.0,ket);
        newstate=auxstate->post_selection(auxprj);
        delete auxprj;
        delete ket;
    }else{
        newstate=auxstate->clone();
    }

    // Free memory
    delete auxstate;

    //Return new state
    return newstate;
}


//----------------------------------------
//
//  Reassign the packet numbers
//
//----------------------------------------
state *state::convert(veci cnv, qocircuit *qoc){
//  veci cnv;         // Index of packets to be reassigned
//  qocircuit *qoc;   // Circuit to which the packets are referred
//  Variables
    int    ch;        // Channel
    int    m;         // Mode
    int    s;         // Packet
    int   *occ;       // Occupation
    state *aux;       // Auxiliary state
    // Auxiliary index
    int    i;         // Aux index
    int    j;         // Aux index
    int    k;         // Aux index


    aux=new state(nlevel,maxket);
    for(i=0;i<nket;i++){
        occ= new int[nlevel]();
        for(j=0;j<nlevel;j++){
            ch=qoc->idx[j].ch;
            m=qoc->idx[j].m;
            if(qoc->idx[j].s<cnv.size()) s=cnv[qoc->idx[j].s];
            else s=qoc->idx[j].s;
            k=qoc->i_idx[ch][m][s];
            occ[k]=occ[k]+ket[i][j];
        }
        aux->add_term(ampl[i],occ);
        delete occ;
    }

    for(j=0;j<nlevel;j++){
        aux->vis[j]=vis[j];
    }

    return aux;
}


//--------------------------------------------
//
// Prints a state in numerical form
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(){


    prnt_state(DEFFORMAT, 0, false, nullptr);
}


//--------------------------------------------
//
// Prints a state in numerical form
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(int column){
//  int column;          // Is the state printed in line or with kets in columns? 0=No/1=Yes


    prnt_state(DEFFORMAT, column, false, nullptr);
}


//--------------------------------------------
//
// Prints a state in human readable form
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(int format, int column, qocircuit *qoc){
//  int column          // Do you want the state print in line or with kets in columns? 0=No/1=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)


    prnt_state(format, column, false,  qoc);
}


//--------------------------------------------
//
// General method to print a state
// ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_state(int format, int column, bool loss, qocircuit *qoc){
//  int format;         // Notation in which the state is printed
//  int column;         // Do you want the state print in line or with kets in columns? 0=No/1=Yes
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)


    if(column==0){
        prnt_in_rows(format, loss,qoc);
    }else{
        prnt_in_cols(format, loss,qoc);
    }

}


//--------------------------------------------
//
// General method to print a state in rows as
// a summation( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_in_rows(int format, bool loss, qocircuit *qoc){
//  int format;         // Notation in which the state is printed
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  Variables
    int  firstline;     // Is this the first term? 1=Yes/0=No
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        if (abs(ampl[i])>DEFTHOLDPRNT){
            if(firstline==0) cout << " + ";

            firstline=0;

            cout << ampl[i];
            cout << " * ";
            prnt_ket(i,format,loss,qoc);
        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//--------------------------------------------
//
// General method to print a state in columns as
// a list ( various kets and amplitudes )
//
//--------------------------------------------
void  state::prnt_in_cols(int format, bool loss, qocircuit *qoc){
//  int format;         // Notation in which the state is printed
//  bool loss;          // Print loss channels in different color? False=No/True=Yes
//  qocircuit *qoc;     // Circuit to which the ket is referred. (Necessary to print in human readable form)
//  Variables
    int  firstline;     // Is this the first term? 1=Yes/0=No
//  Auxiliary index
    int  i;             // Aux index


    firstline=1;
    for(i=0;i<nket;i++){
        if (abs(ampl[i])>DEFTHOLDPRNT){
            prnt_ket(i,format,loss,qoc);


            firstline=0;
            cout << ": ";

            if(real(ampl[i])>=0){
                if(imag(ampl[i])>=0) cout << " " << fixed <<  setprecision(8) << abs(real(ampl[i])) << " + " << abs(imag(ampl[i])) << " j" << endl;
                else                 cout << " " << fixed <<  setprecision(8) << abs(real(ampl[i])) << " - " << abs(imag(ampl[i])) << " j" << endl;
            }else{
                if(imag(ampl[i])>=0) cout << "-" << fixed <<  setprecision(8) << abs(real(ampl[i])) << " + " << abs(imag(ampl[i])) << " j" << endl;
                else                 cout << "-" << fixed <<  setprecision(8) << abs(real(ampl[i])) << " - " << abs(imag(ampl[i])) << " j" << endl;
            }

        }
    }
    if(firstline==1) cout << "| empty >";
    cout << endl;
}


//---------------------------------------------------------------------
//
//  Non-ideal bell emitter with a phase e^(-i phi) in the second term.
//  Auxiliary method.
//  Polarization encoding.
//
//----------------------------------------------------------------------
int state::Bell_Pol(int i_ch0, int i_ch1, veci i_t, char kind, double phi, qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  char   kind;        // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  double phi;         // Non-ideal emission phase
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    double t2;          // Wavepacket 2
    double t3;          // Wavepacket 3
    complex<double> A1; // Amplitude of the first term of the bell state.
    complex<double> A2; // Amplitude of the first term of the bell state.
    mati   Bell1;       // First term of the Bell state
    mati   Bell2;       // Second term of the bell state
    state *auxemitter;  // Auxiliary emitter variable


    // Reserve memory
    Bell1.resize(4,nl);
    Bell2.resize(4,nl);

    // Make variables easily accessible
    t0=i_t(0);
    t1=i_t(1);
    t2=i_t(2);
    t3=i_t(3);

    // Diferent initialization depending of the kind of bell state selected
    switch(kind){
        case '+':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     V,
                         t1,     t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=exp(jm*phi)/sqrt(2.0);

        break;
        case '-':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     V,
                         t1,     t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=-exp(jm*phi)/sqrt(2.0);
        break;
        case 'p':
            Bell1 << i_ch0, i_ch1,
                         H,     V,
                         t0,     t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     H,
                         t1,     t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=exp(jm*phi)/sqrt(2.0);
        break;
        case 'm':
            Bell1 << i_ch0, i_ch1,
                         H,     V,
                         t0,    t2,
                         1,     1;
            Bell2 << i_ch0, i_ch1,
                         V,     H,
                         t1,    t3,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=-exp(jm*phi)/sqrt(2.0);
        break;
        default:
            cout << "Bell emitter error (Polarization): Undefined state" << endl;
            exit(0);
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        if (add_term( A1,Bell1,qoc) <0 ) return -1;
        if (add_term( A2,Bell2,qoc) <0 ) return -1;
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        if ( auxemitter->add_term( A1,Bell1,qoc) <0 ) return -1;
        if ( auxemitter->add_term( A2,Bell2,qoc) <0 ) return -1;
        if ( this->dproduct(auxemitter) <0 ) return -1;
        delete auxemitter;
    }
    return 0;
}


//---------------------------------------------------------------------
//
//  Random noise emitter
//  Auxiliary method.
//
//----------------------------------------------------------------------
int state::Rand_Pol(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    int    r;           // Random integer number between 0 and 3 (included)
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    double t2;          // Wavepacket 2
    double t3;          // Wavepacket 3
    mati   Rand;        // Term of the random state
    state *auxemitter;  // Auxiliary emitter variable


    // Reserve memory
    Rand.resize(4,nl);

    //Make variables easily accessible
    t0=i_t(0);
    t1=i_t(1);
    t2=i_t(2);
    t3=i_t(3);

    // Diferent initialization depending on an aleatory election
    r=(int)(floor(4.0*urand()));
    if(r==4) r=3; // Just in the improbable case we hit the limit.
    switch(r){
        case 0:
            Rand << i_ch0, i_ch1,
                         H,     H,
                         t0,     t2,
                         1,     1;

        break;
        case 1:
            Rand << i_ch0, i_ch1,
                         H,     V,
                         t0,     t3,
                         1,     1;

        break;
        case 2:
            Rand << i_ch0, i_ch1,
                         V,     H,
                         t1,     t2,
                         1,     1;

        break;
        case 3:
            Rand << i_ch0, i_ch1,
                         V,     V,
                         t1,   t3,
                         1,     1;

        break;
        default:
            cout << "Random emitter error: Undefined state" << endl;
            exit(0);
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        if( add_term( 1.0,Rand,qoc) < 0 ) return -1;
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        if( auxemitter->add_term(1.0,Rand,qoc) < 0 ) return -1;
        if( this->dproduct(auxemitter) < 0 ) return -1;
        delete auxemitter;
    }
    return 0;
}


//---------------------------------------------------------------------
//
//  Correlated emitter
//  Auxiliary method.
//
//----------------------------------------------------------------------
int state::Corr_Pol(int i_ch0,int i_ch1, veci i_t, qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
// Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    int    r;           // Random integer number between 0 and 1 (included)
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    double t2;          // Wavepacket 2
    double t3;          // Wavepacket 3
    mati   Corr;        // Term of the correlated state
    state *auxemitter;  // Auxiliary emitter variable


    Corr.resize(4,nl);

    t0=i_t(0);
    t1=i_t(1);
    t2=i_t(2);
    t3=i_t(3);

    // Diferent initialization depending on an aleatory election
    r=(int)(floor(2.0*urand()));
    if(r==2) r=1; // Just in the improbable case we hit the limit.
    switch(r){
        case 0:
            Corr << i_ch0, i_ch1,
                         H,     H,
                        t0,    t2,
                         1,     1;
        break;
        case 1:
            Corr << i_ch0, i_ch1,
                         V,     V,
                        t1,    t3,
                         1,     1;
        break;
        default:
            cout << "Correlated emitter error: Undefined state" << endl;
            exit(0);
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        if( add_term( 1.0,Corr,qoc) < 0 ) return -1;
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        if( auxemitter->add_term( 1.0,Corr,qoc) < 0 ) return -1;
        if( this->dproduct(auxemitter) < 0 ) return -1;
        delete auxemitter;
    }

    return 0;
}


//----------------------------------------------------------------------
//
//  Single pair photon emission in a QD
//  Auxiliary method.
//
//----------------------------------------------------------------------
int state:: QDPair(int i_ch0,int i_ch1, veci i_t,double dt, double k, double S, double tss, double thv, qocircuit *qoc){
//  int    i_ch0;   // Channel 0
//  int    i_ch1;   // Channel 1
//  veci   i_t;     // Vector with the wavepacket numbers to be assigned to each photon.
//  double dt;      // Time between the emission of the two photons.
//  double k;       // fraction of the emitted photon pairs which originate both from a XX-X cascade in presence of background light or multiphoton emission.
//  double S;       // S/hbar FSS constant.
//  double tx;      // Exciton lifetime
//  double tss;     // Time of spin-scattering tss
//  double thv;     // Time of cross-dephasing thv
//  qocircuit *qoc; // Circuit to which the emitter is referred. (Necessary to interpret the defined states)


    if(urand()<(k*exp(-(dt/tss)))){
        if(urand()<exp(-(dt/thv))){
            return Bell_Pol(i_ch0,i_ch1,i_t,'+',S*dt,qoc);
        }else{
            return Corr_Pol(i_ch0,i_ch1,i_t,qoc);
        }
    }else{
        return Rand_Pol(i_ch0,i_ch1,i_t,qoc);
    }

}


//----------------------------------------------------------------------
//
//  Basso Basset Quantum Dot emitter.
//  Simplified version.
//
//----------------------------------------------------------------------
int state:: QD(mati ch, double k, double S, double tx, double tss, double thv, qocircuit *qoc){
//  mati   ch;      // Configuration matrix of the QD. First row channels. Second and third wave packet numbers of the H and V photons respectively.
//  double k;       // Fraction of the emitted photon pairs which originate both from a XX-X cascade in presence of background light or multi-photon emission.
//  double S;       // S/hbar FSS constant.
//  double tx;      // Exciton lifetime
//  double tss;     // Time of spin-scattering tss
//  double thv;     // Time of cross-dephasing thv
//  qocircuit *qoc; // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Variables
    int     ch0;    // Channel 0
    int     ch1;    // Channel 1
    veci    vt;     // Wavepacket numbers vector.
    double  dt;     // Emission time t/tx.
// Note: tx is our code time unit. It is the exciton radiative lifetime.
// Auxiliary index
    int      i;     // Aux index


    vt.resize(4);
    for(i=0;i<ch.cols();i=i+2){
        // Update channel
        ch0  =ch(0,i);
        ch1  =ch(0,i+1);
        vt(0)=ch(1,i);
        vt(1)=ch(2,i);
        vt(2)=ch(1,i+1);
        vt(3)=ch(2,i+1);

        // Emit pair
        dt=tx*expi(urand());
        if ( QDPair(ch0,ch1,vt,dt,k,S,tss,thv,qoc) < 0 ) return -1;
    }
    return 0;
}


//---------------------------------------------------------------------
//
//  Non-ideal bell emitter with a phase e^(-i phi) in the second term.
//  Auxiliary method.
//  Path encoding.
//
//----------------------------------------------------------------------
int state::Bell_Path(int i_ch0, int i_ch1, veci i_t, char kind, double phi,  qocircuit *qoc){
//  int    i_ch0;       // Channel 0
//  int    i_ch1;       // Channel 1
//  veci   i_t;         // Vector with the wavepacket numbers to be assigned to each photon.
//  char   kind;        // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  double phi;         // Non-ideal emission phase
//  qocircuit *qoc;     // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
// Constants
    const  int nl=2;    // Number of level to define for the emitter
//  Variables
    double t0;          // Wavepacket 0
    double t1;          // Wavepacket 1
    complex<double> A1; // Amplitude of the first term of the bell state.
    complex<double> A2; // Amplitude of the first term of the bell state.
    mati   Bell1;       // First term of the Bell state
    mati   Bell2;       // Second term of the bell state
    state *auxemitter;  // Auxiliary emitter variable


    // Reserve memory
    Bell1.resize(4,nl);
    Bell2.resize(4,nl);

    // Make variables easily accessible
    t0=i_t(0);
    t1=i_t(1);
    // Diferent initialization depending of the kind of bell state selected
    switch(kind){
        case '+':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t1,
                         0,     0;
            Bell2 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t1,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=exp(jm*phi)/sqrt(2.0);

        break;
        case '-':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t1,
                         0,     0;
            Bell2 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t1,
                         1,     1;
            A1=1.0/sqrt(2);
            A2=-exp(jm*phi)/sqrt(2.0);
        break;
        case 'p':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t1,
                         0,     1;
            Bell2 << i_ch0, i_ch1,
                         H,     H,
                         t0,     t1,
                         1,     0;
            A1=1.0/sqrt(2);
            A2=exp(jm*phi)/sqrt(2.0);
        break;
        case 'm':
            Bell1 << i_ch0, i_ch1,
                         H,     H,
                         t0,    t1,
                         0,     1;
            Bell2 << i_ch0, i_ch1,
                         H,     H,
                         t0,    t1,
                         1,     0;
            A1=1.0/sqrt(2);
            A2=-exp(jm*phi)/sqrt(2.0);
        break;
        default:
            cout << "Bell emitter error (Path): Undefined state" << endl;
            exit(0);
        break;

    }


    if(nket==0){
        // If the state is not initialized just add the terms of the Bell state
        if ( add_term( A1,Bell1,qoc) < 0 ) return -1;
        if ( add_term( A2,Bell2,qoc) < 0 ) return -1;
    }else{
        // If the state is initialized we have to do the direct product with
        // the already existent term. We are assuming here the pre-existing
        // initialization is for different channels.
        auxemitter=new state(nlevel,maxket);
        if ( auxemitter->add_term( A1,Bell1,qoc) < 0 ) return -1;
        if ( auxemitter->add_term( A2,Bell2,qoc) < 0 ) return -1;
        if ( this->dproduct(auxemitter) < 0 ) return -1;
        delete auxemitter;
    }

    return 0;
}


//---------------------------------------------------------------------
//
//  Non-ideal bell emitter with a phase e^(-i phi) in the second term.
//  Path encoding.
//
//----------------------------------------------------------------------
int state:: Bell(mati ch, char kind, double phi, qocircuit *qoc){
//  mati   ch;      // Configuration matrix of the Bell emitter. First row channels, second row wave packet numbers.
//  char   kind;    // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  double phi;     // Non-ideal emission phase
//  qocircuit *qoc; // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Variables
    int     ch0;    // Channel 0
    int     ch1;    // Channel 1
    veci    vt;     // Wavepacket numbers vector.


    vt.resize(2);

    // Update channel
    ch0  =ch(0,0);
    ch1  =ch(0,1);
    vt(0)=ch(1,0);
    vt(1)=ch(1,1);

    // Emit pair
    return Bell_Path(ch0,ch1,vt,kind,phi,qoc);

}


//---------------------------------------------------------------------
//
//  Non-ideal bell emitter with a phase e^(-i phi) in the second term.
//  Polarization encoding.
//
//----------------------------------------------------------------------
int state:: BellP(mati ch, char kind, double phi, qocircuit *qoc){
//  mati   ch;      // Configuration matrix of the Bell emitter. First row channels. Second and third wave packet numbers of the H and V photons respectively.
//  char   kind;    // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  double phi;     // Non-ideal emission phase
//  qocircuit *qoc; // Circuit to which the emitter is referred. (Necessary to interpret the defined states)
//  Variables
    int     ch0;    // Channel 0
    int     ch1;    // Channel 1
    veci    vt;     // Wavepacket numbers vector.


    vt.resize(4);

    // Update channel
    ch0  =ch(0,0);
    ch1  =ch(0,1);
    vt(0)=ch(1,0);
    vt(1)=ch(1,0);
    vt(2)=ch(1,1);
    vt(3)=ch(1,1);

    // Emit pair
    return Bell_Pol(ch0,ch1,vt,kind,phi,qoc);

}

//----------------------------------------------------------------------
//
//  Direct product of states defined in non coincident channels
//  (Not to be used in general therefore private)
//
//----------------------------------------------------------------------
int state::dproduct(state *rhs){
//  state *rhs      // State in the right hand side to do the d-product
//  Variables
    int    index;   // Ket list position where a new term of the dproduct.
    int   *occ;     // Occupation
    state *aux;     // Auxiliary state
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index


    // Initialize
    aux=new state(this->nlevel,this->maxket);
    aux=this->clone();
    this->clear();
    occ= new int[nlevel]();

    // Dproduct
    for(i=0;i<aux->nket;i++){
        for(j=0;j<rhs->nket;j++){
            for(k=0;k<aux->nlevel;k++){
                occ[k]=aux->ket[i][k]+rhs->ket[j][k];
            }
            index=add_term(aux->ampl[i]*rhs->ampl[j],occ);
            if(index<0){
                cout << "Dproduct: Warning! Memory exceeded. Operation cancelled." << endl;
                // Free memory
                delete[] occ;
                delete aux;
                return -1;
            }
        }
    }

    // Free memory
    delete[] occ;
    delete aux;
    return 0;
}


//---------------------------------------------------------------------------
//
//  Encode state into a qubit representation using path encoding
//
//----------------------------------------------------------------------------
state *state::encode(mati qdef,qocircuit *qoc){
//  mati       qdef; // Integer matrix with the qubit to channel definitions
//  qocircuit *qoc;  // Quatum optical circuit to which this state is referred
//  Variables
    bool   valid;    // Is the ket valid for codification true=Yes/false=No
    bool   print;    // A warning has been printed. True=Yes/False=No
    int    val0;     // value of channel 0
    int    val1;     // value of channel 1
    int    qval;     // qubit equivalent value of channels 0 and 1.
    int    idx;      // Index of the stored state.
    int    nvalid;   // Number of encoded kets
    int   *values;   // Vector of the values of each qubit
    state *qstate;   // qubit encoded state
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index
    int    m;        // Aux index
    int    n;        // Aux index


    // Reserve memory
    qstate=new state(qdef.cols(),maxket);

    // Check encoding conditions
    if((qoc->nm>1)||(qoc->ns>1)){
        cout << "Encode error: The circuit number of modes nm and packets ns can only be one for encoding." << endl;
        return qstate;
    }

    // Encode each ket
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
        // If the resulting state is valid store it
        if(valid==true){
                idx=qstate->add_term(ampl[i],values);
                if((idx!=nvalid)&&(print==false)){
                    cout << "Encode: Warning! encoding leads to collision. Invalid result" << endl;
                    print=true;
                }
                nvalid=nvalid+1;
        }
        delete values;
    }

    // Return encoded state
    return qstate;
}


//---------------------------------------------------------------------------
//
//  Decode the state from a qubit representation into a photon representation in
//  path encoding. State version.
//
//----------------------------------------------------------------------------
state *state::decode(mati qdef,state *ancilla,qocircuit *qoc){
// mati       qdef; // Integer matrix with the qubit to channel definitions
// state *ancilla;  // We need an ancilla state to inform the decoding process of the auxiliary non-qubit channels values
// qocircuit *qoc;  // Quantum optical circuit to which this state is referred
// Variables
    int    val0;    // value of channel 0
    int    val1;    // value of channel 1
    int   *occ;     // Occupation vector
    state *phstate; // Photonic state
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index
    int    l;       // Aux index
    int    m;       // Aux index
    int    n;       // Aux index


    // Reserve and initialize memory
    phstate=new state(ancilla->nlevel,maxket);

    // Check decoding conditions
    if((qoc->nm>1)||(qoc->ns>1)){
        cout << "Decode error: The circuit number of modes nm and packets ns can only be one for encoding." << endl;
        return phstate;
    }

    for(j=0;j<phstate->nlevel;j++) phstate->vis[j]=ancilla->vis[j];

    // Decode each ket
    for(i=0;i<nket;i++){
        occ=new int[phstate->nlevel]();
        for(j=0;j<phstate->nlevel;j++) occ[j]=ancilla->ket[0][j];
        for(j=0;j<qdef.cols();j++){
            // Find the channels
            m=qoc->i_idx[qdef(0,j)][0][0];
            n=qoc->i_idx[qdef(1,j)][0][0];

            // qdef is given over circuit definition.
            // Levels number may change after post-selection
            k=0;
            l=0;
            while((phstate->vis[k]!=m)&&(k<phstate->nlevel)) k=k+1; // Not efficient but small search
            while((phstate->vis[l]!=n)&&(l<phstate->nlevel)) l=l+1; // Not efficient but small search


            // Decode qubit values
            if(ket[i][j]==0){
                val0=0;
                val1=1;
            }else{
                val0=1;
                val1=0;
            }
            occ[k]=val0;
            occ[l]=val1;
        }
        // Store the decoded ket into a new state
        phstate->add_term(ampl[i],occ);
        delete occ;
    }

    // Return result
    return phstate;
}


//---------------------------------------------------------------------------
//
//  Decode the state from a qubit representation into a photon representation in
//  path encoding. Vector version.
//
//----------------------------------------------------------------------------
state *state::decode(mati qdef,veci ancilla,qocircuit *qoc){
//  mati       qdef;     // Integer matrix with the qubit to channel definitions
//  veci ancilla;        // Vector defining the values of the ancilla channels in order ( from smaller to larger channel number )
//  qocircuit *qoc;      // Quantum optical circuit to which this state is referred
//  Variables
    int   *isquch;       // Is a qubit channel? For every channel 0=Ancilla/1= Qubit channel
    mati   def_state;    // State definition
    state *anzstate;     // Ansatz state. State with the ansatz channel values defined.
    state *decoded;      // Decoded state.
//  Auxiliary index
    int    i;            // Aux index
    int    k;            // Aux index


    // Determine which channel is ansatz
    isquch=new int[qoc->nch]();
    for(i=0;i<qdef.cols();i++){
        isquch[qdef(0,i)]=1;
        isquch[qdef(1,i)]=1;
    }

    // Define the state
    k=0;
    def_state.resize(2,qoc->nch);
    for(i=0;i<qoc->nch;i++){
        if(isquch[i]==1){
            // Is not ancilla
            def_state(0,i)=i;
            def_state(1,i)=0;
        }else{
            // Is ancilla
            def_state(0,i)=i;
            def_state(1,i)=ancilla[k];
            k=k+1;
        }
    }

    // Create the ansatz state
    anzstate=new state(qoc->nlevel,1);
    anzstate->add_term(1.0,def_state,qoc);

    // Decode
    decoded=decode(qdef,anzstate,qoc);

    // Free memory.
    delete isquch;
    delete anzstate;

    // Return state
    return decoded;
}


//---------------------------------------------------------------------------
//
//  Encode state into a qubit representation using polarization encoding
//
//----------------------------------------------------------------------------
state *state::pol_encode(veci qdef,qocircuit *qoc){
//  mati       qdef; // Integer matrix with the qubit to channel definitions
//  qocircuit *qoc;  // Quatum optical circuit to which this state is referred
//  Variables
    bool   valid;    // Is the ket valid for codification true=Yes/false=No
    bool   print;    // A warning has been printed. True=Yes/False=No
    int    valH;     // value of channel 0
    int    valV;     // value of channel 1
    int    qval;     // qubit equivalent value of channels 0 and 1.
    int    idx;      // Index of the stored state.
    int    nvalid;   // Number of encoded kets
    int   *values;   // Vector of the values of each qubit
    state *qstate;   // qubit encoded state
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index
    int    m;        // Aux index
    int    n;        // Aux index


    // Reserve memory
    qstate=new state(qdef.size(),maxket);

    // Check encoding conditions
    if((qoc->nm!=2)||(qoc->ns>1)){
        cout << "Pol_encode error: The circuit number of modes nm and packets ns should be nm=2 and ns=1 for encoding." << endl;
        return qstate;
    }

    // Encode each ket
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
        // If the resulting state is valid store it
        if(valid==true){
                idx=qstate->add_term(ampl[i],values);
                if((idx!=nvalid)&&(print==false)){
                    cout << "Pol_encode: Warning! encoding leads to collision. Invalid result" << endl;
                    print=true;
                }
                nvalid=nvalid+1;
        }
        delete values;
    }

    // Return encoded state
    return qstate;
}


//---------------------------------------------------------------------------
//
//  Decode the state from a qubit representation into a photon representation in
//  polarization encoding. State version.
//
//----------------------------------------------------------------------------
state *state::pol_decode(veci qdef,state *ancilla,qocircuit *qoc){
// mati       qdef; // Integer matrix with the qubit to channel definitions
// state *ancilla;  // We need an ancilla state to inform the decoding process of the auxiliary non-qubit channels values
// qocircuit *qoc;  // Quatum optical circuit to which this state is referred
// Variables
    int    valH;    // value of channel 0
    int    valV;    // value of channel 1
    int   *occ;     // Occupation vector
    state *phstate; // Photonic state
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index
    int    k;       // Aux index
    int    l;       // Aux index
    int    m;       // Aux index
    int    n;       // Aux index


    // Reserve and initialize memory
    phstate=new state(ancilla->nlevel,maxket);

    // Check decoding conditions
    if((qoc->nm!=2)||(qoc->ns>1)){
        cout << "Pol_decode error: The circuit number of modes nm and packets ns should be nm=2 and ns=1 for encoding." << endl;
        return phstate;
    }

    for(j=0;j<phstate->nlevel;j++) phstate->vis[j]=ancilla->vis[j];

    // Decode each ket
    for(i=0;i<nket;i++){
        occ=new int[phstate->nlevel]();
        for(j=0;j<phstate->nlevel;j++) occ[j]=ancilla->ket[0][j];
        for(j=0;j<qdef.size();j++){
            // Find the channels
            m=qoc->i_idx[qdef(j)][H][0];
            n=qoc->i_idx[qdef(j)][V][0];

            // qdef is given over circuit definition.
            // Levels number may change after post-selection
            k=0;
            l=0;
            while((phstate->vis[k]!=m)&&(k<phstate->nlevel)) k=k+1; // Not efficient but small search
            while((phstate->vis[l]!=n)&&(l<phstate->nlevel)) l=l+1; // Not efficient but small search


            // Decode qubit values
            if(ket[i][j]==0){
                valH=1;
                valV=0;
            }else{
                valH=0;
                valV=1;
            }
            occ[k]=valH;
            occ[l]=valV;
        }
        // Store the decoded ket into a new state
        phstate->add_term(ampl[i],occ);
        delete occ;
    }

    // Return result
    return phstate;
}


//---------------------------------------------------------------------------
//
//  Decode the state from a qubit representation into a photon representation in
//  polarization encoding. Vector version.
//
//----------------------------------------------------------------------------
state *state::pol_decode(veci qdef,mati ancilla,qocircuit *qoc){
//  veci       qdef;     // Integer vector with the qubit to channel definitions
//  veci ancilla;        // Vector defining the values of the ancilla channels in order ( from smaller to larger channel number )
//  qocircuit *qoc;      // Quantum optical circuit to which this state is referred
//  Variables
    int   *isquch;       // Is a qubit channel? For every channel 0=Ancilla/1= Qubit channel
    mati   def_state;    // State definition
    state *anzstate;     // Ansatz state. State with the ansatz channel values defined.
    state *decoded;      // Decoded state.
//  Auxiliary index
    int    i;            // Aux index
    int    j;            // Aux index
    int    k;            // Aux index
    int    l;            // Aux index


    // Determine which channel is ansatz
    isquch=new int[qoc->nch]();
    for(i=0;i<qdef.size();i++){
        isquch[qdef(i)]=1;
    }

    // Define the state
    k=0;
    l=0;
    def_state.resize(3,qoc->nch*qoc->nm);
    for(i=0;i<qoc->nch;i++){
        for(j=0;j<qoc->nm;j++){
            if(isquch[i]==1){
                // Is not ancilla
                def_state(0,l)=i;
                def_state(1,l)=j;
                def_state(2,l)=0;
            }else{
                // Is ancilla
                def_state(0,l)=i;
                def_state(1,l)=j;
                def_state(2,l)=ancilla(j,k);
                k=k+1;
            }
            l=l+1;
        }
    }

    // Create the ansatz state
    anzstate=new state(qoc->nlevel,1);
    anzstate->add_term(1.0,def_state,qoc);

    // Decode
    decoded=decode(qdef,anzstate,qoc);

    // Free memory.
    delete isquch;
    delete anzstate;

    // Return state
    return decoded;
}


//-------------------------------------------------------
//
//  Creates a projector. The maximum number of kets is set by default.
//
//-------------------------------------------------------
projector::projector(int i_level):state(i_level){
//  int i_level     // Number of levels to describe the projector


    create_projector(i_level,DEFSTATEDIM);
}


//-------------------------------------------------------
//
// Creates a projector specifying the maximum number of kets.
//
//-------------------------------------------------------
projector::projector(int i_level, int i_maxket):state(i_level, i_maxket){
//  int i_level     // Number of levels to describe the projector
//  int i_maxket    // Maximum number of different kets in the summation


    create_projector(i_level,i_maxket);
}


//-----------------------------------------------------
//
//  Creates a projector specifying the maximum number of kets
//  and the vector of equivalence between state levels and
//  circuit defined levels. Intended for internal use.
//
//-----------------------------------------------------
projector::projector(int i_level, int i_maxket, int *i_vis):state(i_level, i_maxket, i_vis){
//  int i_level      // Number of levels to describe the projector
//  int  i_maxket    // Maximum number of different kets in the list
//  int *i_vis       // Vector of equivalence between state levels and circuit defined levels


    create_projector(i_level,i_maxket);
}


//-------------------------------------------------------
//
//  Auxiliary method to create a projector.
//  Same that to create a state but non-defined
//  occupation levels are not assumed to be 0 occupied.
//  They are initialized with a negative value instead.
//  This way these levels are later on ignored/not
// considered in the post-selection operation
//
//-------------------------------------------------------
void projector::create_projector(int i_level, int i_maxket){
//  int i_level     // Number of levels to describe the projector
//  int i_maxket    // Maximum number of different kets in the summation
//  Auxiliary index
    int i;          // Aux index
    int j;          // Aux index


    for(i=0;i<maxket;i++){
        for(j=0;j<nlevel;j++){
            ket[i][j]=-1;
        }
    }
}
