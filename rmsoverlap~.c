/**
 * header file of PD
 * => describes functions, 'classes' and structures used by PD
 */
#include "m_pd.h"
#include <math.h>
#include <string.h>
#include <stdbool.h>
/** CLASS DEFINITION */
static t_class *rmsoverlap_tilde_class;

#define MAXOVERLAP 32
#define INITVSTAKEN 64

/**
 * This is the object's dataspace
 * [1] The first element is necessarily "t_object"
 * [2...] Other data functions as internal variables of the object
 */
/** CLASS STRUCTURE -> variables, memory spaces, inlets and outlets, pointers to various data */
typedef struct _rmsoverlap_tilde {
  t_object x_obj;
  t_float x_f;  
  t_sample x_value;
  t_outlet *x_out;
  t_outlet *x_listout;
  // t_outlet *f_out;
  t_sample *x_win;          /* buffer to define windows (initially, start as 1) */
  t_sample *x_circbuf;      /* circular buffer */  
  t_sample *x_win_overlaps_sum; /* buffer normalization values buffer */  
  t_sample *x_normbuf;
  t_sample *x_sumbuf;
  t_sample *x_layersrms;    /* store here the layers rms val when reach end of window of each layer */
  t_float x_normfact;       /* normalization factor */

  int x_npoints;           /* analysis window size in samples */
  int x_overlaps;          /* overlaps */
  int x_hopsize;  
  int x_size_circbuffer;   /* size in samples of the circular buffer -> in bytes, multiply by sizeof(t_sample) */
  int x_vectorsize;        /* extra buffer for DSP vector size */
  int x_cycles_per_window; /* number of DSP cycles per window */
  int x_dsp_cycles_count;  /* counts from 0 to (x_cycles_per_window - (x_hopsize/x_vectorsize))*/          
  int x_win_init;          /* current sample index where we are writing samples into the circular buffer  */
  char x_winstr[32];
  char x_winalignstr[32];
  char x_normstr[32];
  int x_wintype;
  int x_winaligntype;
  int x_normal_type;    //0: winoverlap_rms
                        //1: winoverlap_mean
                        //2: winoverlap_vals (sample by sample)
                        //3: overlapinverse
                        //4: fixed_multiplier
                        //5: none
  int x_zpadnsamples;
  int *layer_last_cycle;   /* last x_win_init value of each layer */
  int x_fstwinindex;
  int x_padwinsize;
  int x_layer_count;       /* used at the beginning to adjust normalization factor of output according to the number of 
                              analyzed windows: when this reaches the number of overlaps, it will be fixed to that number */
} t_rmsoverlap_tilde;


//AUXILIARY FUNCTIONS
//return next number that is a power of two and greater than arg v
int np2(int x, int v)
{
	float y;
	int out;
	y = pow(2.,ceil(log2f((float)x)));
	if(v > y){
		out = v;
	}else{
		out = (int)y;
	}
  return(out);
}
	
//function to compute windows with zeropaddings and store the result in x->x_win
static void rmsoverlap_tilde_buildwindow(t_rmsoverlap_tilde *x, bool allocwin)
{
  t_sample *win;
  char winstr[32];
  char winalignstr[32];
  int zpadnsamples;
  int npoints;
  int fstwinindex;
  int lstwinindex;
  int padwinsize;


  //voltar daqui!!!
  strcpy(winstr,x->x_winstr);
  strcpy(winalignstr,x->x_winalignstr);
  
  npoints = x->x_npoints;
  zpadnsamples = x->x_zpadnsamples;

  if(allocwin){
    if (!(win = getbytes(sizeof(t_sample) * x->x_size_circbuffer))) 
    {
        pd_error(0, "rmsoverlap~: couldn't allocate internal layers win buffer");
    }
     x->x_win = win;
  }else{
    win = x->x_win;
  } 

  //change wintype according to argument
  if((strcmp(winstr,"rect") == 0) || (strcmp(winstr,"rectangular") == 0) || (strcmp(winstr,"0") == 0)){
    x->x_wintype = 0; //rectangular
    strcpy(x->x_winstr,"rect");
  } else 
  if((strcmp(winstr,"hann") == 0) || (strcmp(winstr,"hanning") == 0) || (strcmp(winstr,"1") == 0)){
     x->x_wintype = 1;
     strcpy(x->x_winstr,"hann");
  } else 
  if((strcmp(winstr,"tri") == 0) || (strcmp(winstr,"triangular") == 0) || (strcmp(winstr,"2") == 0)){
     x->x_wintype = 2;
     strcpy(x->x_winstr,"tri");
  } else 
  if((strcmp(winstr,"hamm") == 0) || (strcmp(winstr,"hamming") == 0) || (strcmp(winstr,"3") == 0)){
     x->x_wintype = 3;
     strcpy(x->x_winstr,"hamm");
  } else 
    if((strcmp(winstr,"black") == 0) || (strcmp(winstr,"blackman") == 0) || (strcmp(winstr,"4") == 0)){
     x->x_wintype = 4;
     strcpy(x->x_winstr,"black");
  } else 
    if((strcmp(winstr,"cos") == 0) || (strcmp(winstr,"cosine") == 0) || (strcmp(winstr,"5") == 0)){
     x->x_wintype = 5;
     strcpy(x->x_winstr,"cos");
  } else {
    x->x_wintype = 0;
    strcpy(x->x_winstr,"rect");
  };  

  // strcpy(winalignstr, x->x_winalignstr);
  //change winalignment according to argument
  if((strcmp(winalignstr,"l") == 0) || (strcmp(winalignstr,"left") == 0) || (strcmp(winalignstr,"0") == 0)){
    x->x_winaligntype = 0; //window aligned to the left side of the analysis block (zeros right-aligned)
    strcpy(x->x_winalignstr,"left");
  } else 
  if((strcmp(winalignstr,"c") == 0) || (strcmp(winalignstr,"center") == 0) || (strcmp(winalignstr,"1") == 0)){
     x->x_winaligntype = 1; //window aligned to the center of analysis block (zeros distributed right and left)
     strcpy(x->x_winalignstr,"center");
  } else 
  if((strcmp(winalignstr,"r") == 0) || (strcmp(winalignstr,"right") == 0) || (strcmp(winalignstr,"2") == 0)){
     x->x_winaligntype = 2; //window aligned to the right side of the analysis block (zeros left-aligned)
     strcpy(x->x_winalignstr,"right");
  } else{
    x->x_winaligntype = 0;
    strcpy(x->x_winalignstr,"left");
  };  

  post("0s %i", zpadnsamples);

  // zeropadding logic
  if ((zpadnsamples > 0)&&(zpadnsamples <= npoints)){
    padwinsize = npoints-zpadnsamples;
    switch(x->x_winaligntype)
    {
      case 0:
        fstwinindex = 0;
          //clean begin with clean window
        for (int i = npoints-zpadnsamples; i < npoints; i++){ //0s at right side
            win[i] = 0.; //right side zero padding
        } 
        break;
      case 1:
        fstwinindex = zpadnsamples/2;
        for (int i = 0; i < fstwinindex; i++){
            win[i] = 0.;           //left side zero padding
            win[npoints-1-i] = 0.; //right side zero padding
        }
        break;
      case 2:
        fstwinindex = zpadnsamples;
        for (int i = 0; i < zpadnsamples; i++){ //0s at left side
            win[i] = 0.; //left side zero padding
        } 
        break;
    }
  }else{
    fstwinindex = 0;
    padwinsize = npoints;
  }

  post("rmsoverlap~: new setup\nwin\t overl\t type\t zpadd\t alig\t norm\n %i \t %i \t %s \t %i \t %s \t %s",  x->x_npoints, x->x_overlaps, x->x_winstr, x->x_zpadnsamples, x->x_winalignstr, x->x_normstr);

  lstwinindex = padwinsize+fstwinindex;
  switch(x->x_wintype)
    {
      case 0: //rectangular
        for (int i = fstwinindex; i < lstwinindex; i++){
          win[i] = 1.; 
        } 
        break;
      case 1: // hanning
        for (int i = fstwinindex; i < lstwinindex; i++){
          win[i] = pow(sin(3.141592653 * (i-fstwinindex) / (padwinsize)),2);
        }
        break;
      case 2: //triangular
        for (int i = fstwinindex, j; i < (padwinsize/2) + fstwinindex; i++){
          j = i+(padwinsize/2);
          win[i] = (float)(i-fstwinindex)/(padwinsize/2);
          win[j] = 1. - (float)(i-fstwinindex)/(padwinsize/2);
        }
        break;      
      case 3: //hamming
        for (int i = fstwinindex; i < lstwinindex; i++){
          win[i] = 0.54-0.46*(cos(6.28318530718*(i-fstwinindex)/padwinsize));
        }
        break;                       
      case 4: //blackman
        for (int i = fstwinindex; i < lstwinindex; i++){
          win[i] = 0.42 - 0.5*(cos(2*3.14159265359*(i-fstwinindex)/padwinsize)) + 0.08*cos((4*3.14159265359*(i-fstwinindex)/padwinsize));
        }
        break;                       
      case 5: //cosine
        for (int i = fstwinindex; i < lstwinindex; i++){
          win[i] = sin(3.14159265359*(i-fstwinindex)/padwinsize);
        }        
        break;  
    } 

  x->x_fstwinindex = fstwinindex;
  x->x_padwinsize = padwinsize;
}

static void rmsoverlap_tilde_normfactor(t_rmsoverlap_tilde *x, int normtype_arg, t_float fixednormfact){

  //0: winoverlap_rms
  //1: winoverlap_mean
  //2: winoverlap_vals (sample by sample)
  //3: overlapinverse
  //4: fixed_multiplier
  //5: none

  // t_sample *win = x->x_win;
  t_sample *win_overlaps_sum;
  t_sample *normbuf;
  // int buf_init_sample;
  int overlaps = x->x_overlaps;
  int npoints = x->x_npoints;
  int win_first_index;
  // int vecsize = x->x_vectorsize;
  // int padwinsize = x->x_padwinsize;
  t_sample normfact;

  normfact = 1.;

  x->x_normal_type = normtype_arg;
  if (!(win_overlaps_sum = getbytes(sizeof(t_sample) * x->x_size_circbuffer))) 
  {
      pd_error(0, "rmsoverlap~: couldn't allocate internal sum buffer");
      // return (0);
  } 
  x->x_win_overlaps_sum = win_overlaps_sum;

  if (!(normbuf = getbytes(sizeof(t_sample) * x->x_size_circbuffer))) 
  {
      pd_error(0, "rmsoverlap~: couldn't allocate internal sum buffer");
      // return (0);
  } 
  x->x_normbuf = normbuf;

  for(int i = 0; i < npoints; i++){
    win_overlaps_sum[i] = 0.;
    normbuf[i] = 0.;
  }
  // Check if any of the overlap layers has ended and store the new RMS value for that layer
  for(int o=0; o < overlaps; o++){
   
    win_first_index= (npoints - (o * x->x_hopsize))%npoints;

    // Apply windowing to the output vector according to the corresponding window segment of each layer
      for(int i = 0; i < npoints; i++){
        int w;
        w = (win_first_index + i)%npoints;
        win_overlaps_sum[i] = win_overlaps_sum[i] + x->x_win[w];
      }
    }
    

  // Important: division by overlaps is calculated in the perform function, which allows it to be dynamic.
  switch(normtype_arg)
  {
    case 5: //none
      for(int i = 0; i < npoints; i++){
        normbuf[i] = 1.;
      }
      break;
    case 4: //fixedval
      for(int i = 0; i < npoints; i++){
        normbuf[i] = fixednormfact;
      }
      break;
    case 3: // divided by overlap
      for(int i = 0; i < npoints; i++){
        normbuf[i] = 1./x->x_overlaps;
      }
      break;      
    case 2: //multiply by inverse of win_overlaps_sum
      for(int i = 0; i < npoints; i++){
        normbuf[i] = 1./win_overlaps_sum[i];
      }
      break;
    case 1:
      //mean of absolute values of the window
      normfact = 0.;
      for(int i = 0; i < npoints; i++){
        normfact += fabs(x->x_win[i]);
      }
      normfact = 1. / ((normfact / npoints) * x->x_overlaps);
      for(int i = 0; i < npoints; i++){
        normbuf[i] = normfact;
      }
      break;

    case 0: //RMS
      normfact = 0.;
      for(int i = 0; i < npoints; i++){
        normfact += (x->x_win[i] * x->x_win[i]);
      }
      normfact = 1. / (sqrt(normfact / npoints) * x->x_overlaps);
      for(int i = 0; i < npoints; i++){
        normbuf[i] = normfact;
      }
      break;
    
  }
  x->x_normfact = normfact;

  
}



/**
 * class constructor function
 */

//static void *rmsoverlap_tilde_new(t_floatarg fnpoints, t_floatarg foverlaps, t_symbol *s)
static void *rmsoverlap_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
  //default values: 1024 2 rect 0. left      
  int npoints = 1024; 
  int overlaps = 2;
  char winstr[32] = "rect"; 
  float padfactor = 0.;
  char winalignstr[32] = "left";
  int normal_type = 0; //win_rms
  int size_circbuffer; 
  int vecsize;
  t_rmsoverlap_tilde *x;
  t_sample *circbuf;
  t_sample *layersrms;
  t_sample *sumbuf;



  x = (t_rmsoverlap_tilde *)pd_new(rmsoverlap_tilde_class);
  x->x_vectorsize = INITVSTAKEN;
  x->x_normal_type = normal_type;
  x->x_layer_count = 0;

  vecsize = x->x_vectorsize;

  switch(argc){
    default:
    case 5: //zero_padding alignment
      atom_string(argv+4, winalignstr, 10);
    case 4: //zero_padding (number of samples padded to 0)
      padfactor = atom_getfloat(argv+3);
    case 3:
      atom_string(argv+2, winstr, 20);
    case 2:
      overlaps = atom_getfloat(argv+1);
    case 1:
      npoints = atom_getfloat(argv);
  }

  strcpy(x->x_winstr, winstr); 
  strcpy(x->x_winalignstr, winalignstr); 
 

  // set default values for npoints and overlaps
  if (npoints < 1){
    npoints = 1024;
  }else{
    npoints = np2(npoints, x->x_vectorsize);
  }; 
  
  if (overlaps < 1){
    overlaps = 1;
  }else if (overlaps >= npoints ){
    overlaps = npoints;
  }else{
    overlaps =  np2(overlaps, 1);
  }


    post("analysis block: %i, overlaps: %i",npoints,overlaps);
    
    //  post("npoints: %i", npoints);
    x->x_npoints = npoints; //np2(npoints, x->x_vectorsize);
    x->x_overlaps = overlaps;
    x->x_zpadnsamples = (int)(padfactor * npoints);

    
  // circular buffer size in number of samples (for bytes multiply by sizeof(t_sample))
  size_circbuffer = (x->x_npoints + x->x_vectorsize);
  x->x_size_circbuffer = size_circbuffer; 
  
  // allocation of circular buffer
  if (!(circbuf = getbytes(sizeof(t_sample) * size_circbuffer))) 
  {
    pd_error(0, "rmsoverlap~: couldn't allocate internal circular buffer");
    return (0);
  }    

  // allocation of buffer for RMS value of each overlap layer 
    if (!(layersrms = getbytes(sizeof(t_sample) * overlaps))) 
    {
        pd_error(0, "rmsoverlap~: couldn't allocate internal layers rms buffer");
        return (0);
    }    

    if (!(sumbuf = getbytes(sizeof(t_sample) * vecsize))) 
    {
        pd_error(0, "rmsoverlap~: couldn't allocate internal sum buffer");
        return (0);
    } 
    x->x_sumbuf = sumbuf;

    x->x_circbuf = circbuf;
    x->x_layersrms = layersrms;
    x->x_hopsize = x->x_npoints/x->x_overlaps;
    x->x_cycles_per_window = x->x_npoints/x->x_vectorsize;
    x->x_dsp_cycles_count = 0;
    x->x_win_init = 0;
    
    //zero circular buffer at the beginning
    for(int i=0; i < size_circbuffer; i++){
      circbuf[i] = 0.;
    }   

    // Initialize RMS for each layer to 0: possibility -> use negative value and, from that, change normalization
    for(int i=0; i < overlaps;i++){
     layersrms[i] = 0.;
    }

    for(int i=0; i < vecsize;i++){
     sumbuf[i] = 0.;
    }

  //auxiliary function to build window
  rmsoverlap_tilde_buildwindow(x, true);
  strcpy(x->x_normstr,"rect");
  rmsoverlap_tilde_normfactor(x, 0, 0.); //default normalization: rms of window-overlaps

  /* create a new signal-outlet */
  x->x_out = outlet_new(&x->x_obj, &s_signal);
  x->x_listout = outlet_new(&x->x_obj, gensym("list"));

  return (x);
}

static void rmsoverlap_tilde_free(t_rmsoverlap_tilde *x)
{
  /* free any ressources associated with the given outlet */
  freebytes(x->x_circbuf, sizeof(t_sample) * x->x_size_circbuffer);
  freebytes(x->x_layersrms, sizeof(t_sample) * x->x_overlaps);
  freebytes(x->x_win, sizeof(t_sample) * x->x_size_circbuffer);
  freebytes(x->x_win_overlaps_sum, sizeof(t_sample) * x->x_size_circbuffer);
  freebytes(x->x_normbuf, sizeof(t_sample) * x->x_size_circbuffer);
  freebytes(x->x_sumbuf, sizeof(t_sample) * x->x_vectorsize);
  outlet_free(x->x_out);
  outlet_free(x->x_listout);
}



//rever
static t_int *rmsoverlap_tilde_perform(t_int *w)
{
    /* [1] pointer to the data space of the object itself */
    t_rmsoverlap_tilde *x = (t_rmsoverlap_tilde *)(w[1]);
    /* [2] pointer to the array of type t_sample that points to the input signal */
    t_sample *in = (t_sample *)(w[2]);
    /* [3] pointer to the array of type t_sample that points to the output signal */
    t_sample *out = (t_sample *)(w[3]);
    /* [4] size of the input vector */
    int vecsize = w[4];
    t_sample *readp = in;
    t_sample *cbuf = x->x_circbuf;
    t_sample *layersrms = x->x_layersrms;
    t_sample *win = x->x_win;
    t_sample *sumbuf;
    t_sample *win_overlaps_sum = x->x_win_overlaps_sum;
    t_sample *normbuf = x->x_normbuf;
    t_float normfact = x->x_normfact;
    int nextwininit;
    int buf_init_sample;
    int overlaps = x->x_overlaps;
    int npoints = x->x_npoints;
    int win_first_index;
    float overlap_norm = 1.;

    

    
    buf_init_sample = x->x_win_init;
    nextwininit = (x->x_win_init + vecsize)%(x->x_npoints);
    sumbuf = x->x_sumbuf;

   

    // //write new cycle of incoming vector to circular buffer
    for(int i=0; i<vecsize; i++){
      int j;
      j = (i + buf_init_sample) % (x->x_npoints); // j: index of the start of writing the current vector in the circular buffer
      
      cbuf[j] = readp[i];
      sumbuf[i] = 0.;
    }

    //post("buf_init_sample: %i \t nextwininit: %i", buf_init_sample, nextwininit);  
    

    // check if any of the overlap layers has ended and store the new RMS value for that layer
    for(int o=0; o < overlaps; o++){

      //if layer will end at the next start point, calculate new value for the layer
      if(nextwininit == ((float)npoints * o)/(float)overlaps){
       // post("o: %i \t || nextwininit: %i || match: %f", o, nextwininit, ((float)npoints * o)/(float)overlaps);  
        t_sample rms = 0.;
        
        if(x->x_layer_count < overlaps){
          x->x_layer_count += 1;
        }

        for(int s=0; s < npoints; s++){
          int j;
          j = (s+nextwininit)%(npoints);
          rms = rms + (cbuf[j] * cbuf[j]);
        }
        rms = sqrt(rms/(float)npoints);
        layersrms[o] = rms;
      }
    
    win_first_index= (npoints - (o * x->x_hopsize) + buf_init_sample)%npoints;

      for(int i = 0; i < vecsize; i++){
        int wi;
        wi = (win_first_index+ i)%npoints;// (i + (x->x_hopsize * (x->x_overlaps - o))) % x->x_npoints;
        //post("o, i, w: %i  %i  %i", o, i, w);
        sumbuf[i] = sumbuf[i] + (win[wi] * layersrms[o]);
      }
    }


  //0: winoverlap_rms
  //1: winoverlap_mean
  //2: winoverlap_vals (sample by sample)
  //3: overlapinverse
  //4: fixed_multiplier
  //5: none

  for(int i = 0; i < vecsize; i++){
    out[i] = sumbuf[i] * normbuf[i];
  }


    x->x_win_init = nextwininit;
    x->x_dsp_cycles_count = x->x_win_init/vecsize;
    return (w+5); // precisar ser N+1 (para N argumentos)
}


/* executado quando liga-se o DSP */
/**
 * Resets the size of a circular buffer and other parameters when DSP is turned on.
 *
 * @param x     A pointer to the rmsoverlap_tilde object.
 * @param sp    An array of signal vectors.
 *              sp[0]->s_vec: Input vector
 *              sp[1]->s_vec: Output vector
 * @remarks     This function is executed when the DSP is turned on for the rmsoverlap~ object.
 *              It checks if the block size has changed and resizes the circular buffer if necessary.
 *              Then it adds the rmsoverlap_tilde_perform function to the DSP chain for audio processing.
 */
static void rmsoverlap_tilde_dsp(t_rmsoverlap_tilde *x, t_signal **sp)
{
  size_t new_circbuffer_size;
  t_int blocksize;

  blocksize = (t_int)sp[0]->s_n;

  // confere se mudou block size
  if (blocksize != (t_int)x->x_vectorsize) // era "<".. troquei por "!="
    {
      new_circbuffer_size = sizeof(t_sample) * (x->x_npoints);
      //pointer to void: 
      //  maneira de especificar ponteiro sem definir tipo do dado armazenado na memória
      //  em um segundo momento, usa-se typecast para fazer o dado ser interpretadode uma maneira específica
      //a função resizebytes é usada para realocar a memória de x_circbuf caso, ao ligar o dsp,
      //detecte-se que o espaço anteriormente alocado não é suficiente para armazenar os dados
      void *xx = resizebytes(x->x_circbuf, // parece ser necessário caso mude o tamanho do vetor...
          x->x_size_circbuffer,
          new_circbuffer_size);
      if (!xx)
      {
          pd_error(0, "rmsoverlap~: out of memory");
          return;
      }
      x->x_circbuf = (t_sample *)xx;
      x->x_vectorsize = blocksize;
      x->x_cycles_per_window = x->x_npoints / x->x_vectorsize; ///algo, anteriormente, para verificar potencias de 2?
    }

    dsp_add(rmsoverlap_tilde_perform, 
            4,
            x,                  //objeto
            sp[0]->s_vec,       //input vector
            sp[1]->s_vec,       //output vector
            (t_int)sp[0]->s_n   //vector_size       
            );
            
}


void rmsoverlap_tilde_window(t_rmsoverlap_tilde *x, t_symbol *sel, int argc, t_atom *argv)
{ 
  char winstr[32];  

  // winstr = x->x_winstr;

  if (argv[0].a_type == A_FLOAT){
    winstr[0] = (int)(argv->a_w.w_float + 48);
  };

  if (argv[0].a_type == A_SYMBOL){
    strcpy(winstr,argv->a_w.w_symbol->s_name);
  };

  strcpy(x->x_winstr, winstr);
  // x->x_winstr = winstr;
  
  rmsoverlap_tilde_buildwindow(x, false);
}

void rmsoverlap_tilde_windowalign(t_rmsoverlap_tilde *x, t_symbol *sel, int argc, t_atom *argv)
{ 
  char winstralign[32];  

  // winstr = x->x_winstr;

  if (argv[0].a_type == A_FLOAT){
    winstralign[0] = (int)(argv->a_w.w_float + 48);
  };

  if (argv[0].a_type == A_SYMBOL){
    strcpy(winstralign,argv->a_w.w_symbol->s_name);
    post("===> %s", argv->a_w.w_symbol->s_name);
  };

  strcpy(x->x_winalignstr, winstralign);
  // x->x_winalignstr = winstralign;
  
  rmsoverlap_tilde_buildwindow(x, false);
  
}


void rmsoverlap_tilde_zeropadding(t_rmsoverlap_tilde *x, t_symbol *sel, int argc, t_atom *argv)
{ 
  float zeropaddval;
  int zeros;

  if (argv[0].a_type == A_FLOAT){
    zeropaddval = (argv->a_w.w_float);
    zeros = (int)(zeropaddval * x->x_npoints);
    x->x_zpadnsamples = zeros;

    rmsoverlap_tilde_buildwindow(x, false);
    // post("rmsoverlap~: new setup\nwindow\t overlaps\t type\t zeropadding\t alignment\n %i \t %i \t %s \t %i \t %s",  x->x_npoints, x->x_overlaps, x->x_winstr, x->x_zpadnsamples, x->x_winalignstr);
  }; 
}


void rmsoverlap_tilde_getwindow(t_rmsoverlap_tilde *x, t_symbol *s)
{
  int npoints;
  t_atom *atombuf;
  t_sample *win;
   
  npoints = x->x_npoints;
  atombuf = (t_atom *)getbytes(sizeof(t_atom)*(npoints));
  win = x->x_win;
  
  for (int n = 0; n < npoints; n++) {
    SETFLOAT(&atombuf[n], win[n]);
  }

  outlet_anything(x->x_listout, gensym("window"), npoints, atombuf);
  
  freebytes(atombuf,sizeof(t_atom)*npoints);
}

void rmsoverlap_tilde_getwindowsum(t_rmsoverlap_tilde *x, t_symbol *s)
{
  int npoints;
  t_atom *atombuf;
  t_sample *win_overlaps_sum;
   
  npoints = x->x_npoints;
  atombuf = (t_atom *)getbytes(sizeof(t_atom)*(npoints));
  win_overlaps_sum = x->x_win_overlaps_sum;
  
  for (int n = 0; n < npoints; n++) {
    SETFLOAT(&atombuf[n], win_overlaps_sum[n]);
  }

  outlet_anything(x->x_listout, gensym("window_sum"), npoints, atombuf);
  
  freebytes(atombuf,sizeof(t_atom)*npoints);
}

void rmsoverlap_tilde_getnormvals(t_rmsoverlap_tilde *x, t_symbol *s)
{
  int npoints;
  t_atom *atombuf;
  t_sample *normbuf;
   
  npoints = x->x_npoints;
  atombuf = (t_atom *)getbytes(sizeof(t_atom)*(npoints));
  normbuf = x->x_normbuf;
  
  for (int n = 0; n < npoints; n++) {
    SETFLOAT(&atombuf[n], normbuf[n]);
  }

  outlet_anything(x->x_listout, gensym("window_norm"), npoints, atombuf);
  
  freebytes(atombuf,sizeof(t_atom)*npoints);
}





void rmsoverlap_tilde_normalize(t_rmsoverlap_tilde *x, t_symbol *sel, int argc, t_atom *argv)
{ 
    char normstr[32];
  
    // post("normalize: selector %s", sel->s_name);

      //0: winoverlap_rms
      //1: winoverlap_mean
      //2: winoverlap_vals (sample by sample)
      //3: overlap_inverse
      //4: fixed_multiplier
      //5: none

      //if number, apply the value as a fixed output multiplier
    	if (argv[0].a_type == A_FLOAT){
        x->x_normfact = argv[0].a_w.w_float;
        rmsoverlap_tilde_normfactor(x, 4, argv[0].a_w.w_float);
	      post("rmsoverlap: output vals will be multiplied by %.8f", argv[0].a_w.w_float);
      } else if (argv[0].a_type == A_SYMBOL){
        strcpy(normstr,argv[0].a_w.w_symbol->s_name);
        if(strcmp(normstr,"none") == 0){
          rmsoverlap_tilde_normfactor(x, 5, 1.);
        }else 
        if(strcmp(normstr,"fixed") == 0){
          rmsoverlap_tilde_normfactor(x, 4, (float)argv[1].a_w.w_float);
        }else 
        if(strcmp(normstr,"overlap") == 0){
          rmsoverlap_tilde_normfactor(x, 3, 0.);
        }else 
        if(strcmp(normstr,"win_vals") == 0){
          rmsoverlap_tilde_normfactor(x, 2, 0.);
        }else
        if(strcmp(normstr,"win_mean") == 0){
          rmsoverlap_tilde_normfactor(x, 1, 0.);
        }else
        if(strcmp(normstr,"win_rms") == 0){
          rmsoverlap_tilde_normfactor(x, 0, 0.);
        }else{
	      post("rmsoverlap~: normalization option not recognized.\nWill use win_rms.");
         rmsoverlap_tilde_normfactor(x, 0, 0.);
      }
      }
      strcpy(x->x_normstr,normstr);
}



void rmsoverlap_tilde_setup(void)
{
    rmsoverlap_tilde_class = 
        class_new(gensym("rmsoverlap~"),   // [1] 1o argumento: nome simbólico da classe
        (t_newmethod)rmsoverlap_tilde_new, // [2] 2o argumento: função de construção
        (t_method)rmsoverlap_tilde_free,   // [3] 3o argumento: função de destrução         
        sizeof(t_rmsoverlap_tilde),        // [4] 4o argumento: Para permitir ao Pd reservar e liberar memória 
                                        //     suficiente para o dataspace estático, o tamanho da etrutura 
                                        //     de dados precisa ser informada como 4o argumento
        CLASS_DEFAULT,                  // [5] 5o argumento: influencia a representação gráfica dos objetos da classe.
                                        //     O valor default é CLASS_DEFAULT ou simplesmente "0".
                                               
        A_GIMME,                     // [6...] Os argumento restantes de um objeto e seu tipo. Até 6 argumentos numéricos
                            //     (A_DEFFLOAT) ou simbólicos (A_DEFSYMBOL) podem ser definidos.
                          //     Se for necessário passar mais de 6 argumentos, A_GIMME pode ser utilizado
                                        //     para passar uma lista de argumentos. 
                                        
                                        //  Aqui não está indicado nenhum argumento!
                                        
                               
        0);                             // A lista de argumentos do objeto é terminada com “0”. 
        
    CLASS_MAINSIGNALIN(rmsoverlap_tilde_class, t_rmsoverlap_tilde, x_f);
    class_addmethod(rmsoverlap_tilde_class, (t_method)rmsoverlap_tilde_dsp,
        gensym("dsp"), A_CANT, 0);

    class_addmethod(rmsoverlap_tilde_class,
        (t_method)rmsoverlap_tilde_normalize, gensym("normalize"),
        A_GIMME, 0);

    class_addmethod(rmsoverlap_tilde_class,
        (t_method)rmsoverlap_tilde_window, gensym("window"),
        A_GIMME, 0);        

    class_addmethod(rmsoverlap_tilde_class,
        (t_method)rmsoverlap_tilde_windowalign, gensym("winalign"),
        A_GIMME, 0);           

    class_addmethod(rmsoverlap_tilde_class,
        (t_method)rmsoverlap_tilde_zeropadding, gensym("zeropadding"),
        A_GIMME, 0);

    class_addmethod(rmsoverlap_tilde_class,
      (t_method)rmsoverlap_tilde_getwindow, gensym("get_window"),
      A_GIMME, 0);      

    class_addmethod(rmsoverlap_tilde_class,
      (t_method)rmsoverlap_tilde_getwindowsum, gensym("get_window_sum"),
      A_GIMME, 0);

     class_addmethod(rmsoverlap_tilde_class,
      (t_method)rmsoverlap_tilde_getnormvals, gensym("get_norm_vals"),
      A_GIMME, 0);    

}

