/*** Code for generation of baseband and modulated passband waveforms for multiple symbols. ***/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unordered_map>
#include <cmath>
#include <vector>
using namespace std;

unordered_map<int, string> umap;
vector<float> arr_even_x;
vector<float> arr_even_y;
vector<float> arr_even_y_imag;
vector<float> arr_odd_x;
vector<float> arr_odd_y;
vector<float> arr_odd_y_imag;
vector<float> com_arr_y;
vector<float> com_arr_y_imag;
double pi = 2 * acos(0.0);

/*** Function Definition for Bits to Symbol Mapping***/

int btsm(vector<int> &arr)
{
    /*** Converting 4 bit to decimal symbol value (0-15) ***/

    int btsm_val = 0;
    int arr_size = arr.size();
    for (int i = 0; i < arr_size; i++)
    {
        btsm_val += pow(2, i) * arr[i];
    }
    return btsm_val;
}

/*** Function Definition for Symbol to Chip Mapping ***/

string stcm(int sym)
{
    /***Using unordered map for storing the chip sequences corresponding to symbol.***/

    umap[0] = "11011001110000110101001000101110";
    umap[1] = "11101101100111000011010100100010";
    umap[2] = "00101110110110011100001101010010";
    umap[3] = "00100010111011011001110000110101";
    umap[4] = "01010010001011101101100111000011";
    umap[5] = "00110101001000101110110110011100";
    umap[6] = "11000011010100100010111011011001";
    umap[7] = "10011100001101010010001011101101";
    umap[8] = "10001100100101100000011101111011";
    umap[9] = "10111000110010010110000001110111";
    umap[10] = "01111011100011001001011000000111";
    umap[11] = "01110111101110001100100101100000";
    umap[12] = "00000111011110111000110010010110";
    umap[13] = "01100000011101111011100011001001";
    umap[14] = "10010110000001110111101110001100";
    umap[15] = "11001001011000000111011110111000";
    return umap[sym];
}

/*** Function and struct definition for splitting the chip sequences into the odd and even numbered chip sequences. ***/

struct e_o_strings
{
    string even;
    string odd;
};

typedef struct e_o_strings Struct;

Struct c_seq_split(string c_map)
{
    Struct s;
    
    /*** Separating into even and odd numbered chip sequences in the form of strings. ***/
    
    for (int i = 0; i < c_map.length(); i++) {
        if (i % 2 == 0)
            s.even.push_back(c_map[i]);
        else
            s.odd.push_back(c_map[i]);
    }
    return s;
}

/*** Function to generate the quadrature component of the OQPSK baseband waveform per data symbol. ***/
/*** Takes the even chip sequence, Chip period and updated count variable and Samples Per Chip as parameters. ***/

int oqpsk_ip(string c_even_map, float &Tc, int &count, int &spc)
{
    /*** The step_size parameter decided by the samplesPerChip parameter for storing the time t values in each chip pulse. ***/
    float step_size = Tc/spc; //Step_size = Chip Period / Samples Per Chip (Time slot of each sample in each chip pulse)
    for (int i = 0; i < c_even_map.length(); i++)
    {
        /***Iterating over the even chip sequence and performing pulse shaping based on individual chip values.***/
        
        if (c_even_map[i] == '1')
        {
        /*** Inside if statement, looping from 0 to 2 * spc, every time 1 or 0 is obtained as the chip value and 
         applying the corresponding half sin pulse shaping. The corresponding x and y values are stored in different 
         vectors which are then used for plotting. count variable tracks the number of chips visited by the loop. ***/
         
        /***The loop is run from 0 to 2 * spc for each chip pulse. ***/

            for (int k = 0; k <= (2 * spc); k++)
            {
            /*** The multiplication by step_size to k is done to get the actual value of time t for pulse shaping. ***/
                arr_even_x.push_back(k * step_size + (count));
                if (abs(sin(pi * k * step_size / (2 * Tc))) < 1e-6) // Inserting zeroes for very small sin values.
                    arr_even_y.push_back(0);
                else
                    arr_even_y.push_back(sin(pi * k * step_size/(2 * Tc)));
            }
        }
        else
        {
            for (int k = 0; k <= (2 * spc);k++)
            {
                arr_even_x.push_back(k * step_size + (count));
                if (abs(-1 * sin(pi * k * step_size / (2 * Tc))) < 1e-6)
                    arr_even_y.push_back(0);
                else
                    arr_even_y.push_back(-1 * sin(pi * k * step_size / (2 * Tc)));
            }
        }
        count += 2 * Tc;  // Moving count to the next pulse.
    }

    /*** Commented code for printing the x,y values eventually to be used for plotting. ***/
    
    /*cout << "       " << endl;
    cout << "ip comp ----" << endl;
    cout << "       " << endl;
    for (int i = 0; i < arr_even_x.size(); i++)
    {
        cout << arr_even_x[i] << ',' << arr_even_y[i] << endl;
    }*/
    return count; //Returning count variable for multiple symbols case to keep correct track of the time on x axis.
}

/*** Function to generate the quadrature component of the OQPSK baseband waveform per data symbol. ***/
/*** Takes the odd chip sequence, Chip period, updated count variable and Samples Per Chip as parameters. ***/

int oqpsk_qp(string c_odd_map, float &Tc, int &count, int &spc)  
{
    /*** Iterating over the odd chip sequence and performing pulse shaping based on individual chip values.* ***/
    
    float step_size = Tc/spc;
    bool flag = true; // boolean variable to insert initial x,y values from 0 to Tc in case of quadrature phase.
    for (int i = 0; i < c_odd_map.length(); i++)
    {
        /*** Inside if statement, looping from spc to 3 * spc, every time 1 or 0 is obtained as the chip value and 
         applying the corresponding half sin pulse shaping. The corresponding x and y values are stored in different 
         vectors which are then used for plotting. count variable tracks the number of chips visited by the loop. ***/

        if (c_odd_map[i] == '1')
        {
            /***The loop is run from spc to 3 * spc for each chip pulse.***/

            for (int k = spc; k <= 3 * spc; k++)
            {
                /*** The logic to add the offset at the beginning of the quadrature phase sequence. ***/
                if (count == 0 && flag)
                {
                    int m = 0;
                    while (m < spc)
                    {
                        arr_odd_x.push_back(m * step_size);
                        arr_odd_y.push_back(0);
                        m++;
                    }
                }
                flag = false; //setting the boolean flag false after the first Tc values on the first count are recorded.
                
                /*** The multiplication by step_size of k is done to get the actual value of time t for pulse shaping. ***/

                arr_odd_x.push_back(k * step_size + count);
                if (abs(sin(pi * (k * step_size - Tc) / (2 * Tc))) < 1e-6) // Inserting zeroes for very small sin values.
                    arr_odd_y.push_back(0);
                else
                    arr_odd_y.push_back(sin(pi * (k * step_size - Tc) / (2 * Tc)));
            }
        }
        else
        {
            for (int k = spc; k <= 3 * spc; k++)
            {
                if (count == 0 && flag)
                {
                    int m = 0;
                    while (m < spc)
                    {
                        arr_odd_x.push_back(m * step_size);
                        arr_odd_y.push_back(0);
                        m++;
                    }
                }
                flag = false;
                arr_odd_x.push_back(k * step_size + count);
                if (abs(-1 * sin(pi * (k * step_size - Tc) / (2 * Tc))) < 1e-6)
                    arr_odd_y.push_back(0);
                else
                    arr_odd_y.push_back(-1 * sin(pi * (k * step_size - Tc) / (2 * Tc)));
            }
        }
        count += 2 * Tc;  // Moving count to the next pulse.
    }

    /*** Commented code for printing the x,y values eventually to be used for plotting. ***/
    cout << "       " << endl;
    cout << "qp comp ----" << endl;
    cout << "       " << endl;
    for (int i = 0; i < arr_odd_x.size(); i++)
    {
        cout << arr_odd_x[i] << ',' << arr_odd_y[i] << endl;
    }
    cout << arr_odd_x.size() << endl;
    return count; //Returning count variable for multiple symbols case to keep correct track of the time on x axis.
}

/*** Function to generate the I(t)*Cos(wct) and I(t)*Sin(wct) terms of OQPSK modulated passband waveform. ***/
/*** Takes even chip sequence, chip period Tc, carrier frequency fc, updated count variable and Samples Per Chip as parameters. ***/

int com_ip_oqpsk(string &c_even_map, float &Tc, long long &fc, int &count, int &spc)
{
    float step_size = Tc/spc;
    float wc = 2 * pi * fc; //Calculation of wc
    for (int i = 0; i < c_even_map.length(); i++)
    {
        if (c_even_map[i] == '1')
        {
            for (int k = 0; k <= (2 * spc); k++)
            {
                /*** Pushing back the value k*step_size + count for x values, to get the correct value of time t.  
                 Multiplying the value of time t with 1e-6 while computing Cos(wct) and Sin(wct) as t is measured in us. ***/
                arr_even_x.push_back(k * step_size + (count));
                if (abs(sin(pi * k * step_size / (2 * Tc))) < 1e-6) {
                    arr_even_y.push_back(0);
                    arr_even_y_imag.push_back(0);
                }
                else {
                    arr_even_y.push_back(sin(pi * k * step_size / (2 * Tc)) * (cos(wc * (k * double(step_size) + (count)) * 1e-6)));
                    arr_even_y_imag.push_back(sin(pi * k * step_size / (2 * Tc)) * (sin(wc * (k * double(step_size) + (count)) * 1e-6)));
                }
            }
        }
        else
        {
            for (int k = 0; k <= (2 * spc); k++)
            {
                arr_even_x.push_back(k * step_size + (count));
                if (abs(-1 * sin(pi * k * step_size / (2 * Tc))) < 1e-6) {
                    arr_even_y.push_back(0);
                    arr_even_y_imag.push_back(0);
                }
                else {
                    arr_even_y.push_back(-1 * sin(pi * k * step_size / (2 * Tc)) * (cos(wc * (k * double(step_size) + (count)) * 1e-6)));
                    arr_even_y_imag.push_back(-1 * sin(pi * k * step_size / (2 * Tc)) * (sin(wc * (k * double(step_size) + (count)) * 1e-6)));
                }
            }
        }
        count += 2 * Tc;
    }
    /*** Commented code for printing the x,y values eventually to be used for plotting. ***/

    /*for (int i=0;i<arr_even_y.size();i++){
        cout << arr_even_x[i] << ',' << arr_even_y[i] << endl;
    }*/
    return count; //Returning count variable for multiple symbols case to keep correct track of the time on x axis.
}

/*** Function to generate the Q(t)*Sin(wct) and Q(t)*Cos(wct) terms of OQPSK modulated passband waveform. ***/
/*** Takes odd chip sequence, chip period Tc, carrier frequency fc, updated count variable and Samples Per Chip as parameters. ***/

int com_qp_oqpsk(string &c_odd_map, float &Tc, long long &fc, int &count, int &spc)
{
    float wc = 2 * pi * fc;
    float step_size = Tc/spc;
    bool flag = true;
    for (int i = 0; i < c_odd_map.length(); i++)
    {
        if (c_odd_map[i] == '1')
        {
            for (int k = spc; k <= 3 * spc; k++)
            {
                if (count == 0 && flag)
                {
                    int m = 0;
                    while (m < spc)
                    {
                        arr_odd_x.push_back(m * step_size);
                        arr_odd_y.push_back(0);
                        arr_odd_y_imag.push_back(0);
                        m++;
                    }
                }
                flag = false;
                /*** Pushing back the value k*step_size + count for x values, to get the correct value of time t.  
                 Multiplying the value of time t with 1e-6 while computing Sin(wct) and Cos (wct) as t is measured in us. ***/
                arr_odd_x.push_back(k * step_size + count);
                if (abs(sin(pi * (k * step_size - Tc) / (2 * Tc))) < 1e-6) {
                    arr_odd_y.push_back(0);
                    arr_odd_y_imag.push_back(0);
                }
                else {
                    arr_odd_y.push_back(sin(pi * (k * step_size - Tc) / (2 * Tc)) * sin(wc * (k * double(step_size) + count) * 1e-6));
                    arr_odd_y_imag.push_back(sin(pi * (k * step_size - Tc) / (2 * Tc)) * cos(wc * (k * double(step_size) + count) * 1e-6));
                }
            }
        }
        else
        {
            for (int k = spc; k <= 3 * spc; k++)
            {
                if (count == 0 && flag)
                {
                    int m = 0;
                    while (m < spc)
                    {
                        arr_odd_x.push_back(m * step_size);
                        arr_odd_y.push_back(0);
                        arr_odd_y_imag.push_back(0);
                        m++;
                    }
                }
                flag = false;
                arr_odd_x.push_back(k * step_size + count);
                if (abs(-1 * sin(pi * (k * step_size - Tc) / (2 * Tc))) < 1e-6){
                    arr_odd_y.push_back(0);
                    arr_odd_y_imag.push_back(0);
                }    
                else {
                    arr_odd_y.push_back(-1 * sin(pi * (k * step_size - Tc) / (2 * Tc)) * sin(wc * (k * double(step_size) + count) * 1e-6));
                    arr_odd_y_imag.push_back(-1 * sin(pi * (k * step_size - Tc) / (2 * Tc)) * cos(wc * (k * double(step_size) + count) * 1e-6));
                }
            }
        }
        count += 2 * Tc;
    }
    /*** Commented code for printing the x,y values eventually to be used for plotting. ***/

    /*for (int i=0;i<arr_odd_x.size();i++){
        cout << arr_odd_x[i] << ',' << arr_odd_y[i] << endl;
    }*/
    return count; //Returning count variable for multiple symbols case to keep correct track of the time on x axis.
}

/*** Function to subtract the Q(t)*Sin(wct) term from the I(t)*Cos(wct) term to get the real part of 
     OQPSK modulated passband waveform. ***/

void real_part_com_ip_qp(vector<float> &arr_even_x, vector<float> &arr_even_y, vector<float> &arr_odd_y)
{
    /*** Values of the I(t)Cos(wct) and Q(t)Sin(wct) computed from the prior functions are subtracted from each other 
     to compute the combined real part of the passband modulated waveform in time domain. ***/

    for (int i = 0; i < arr_even_x.size(); i++)
    {
        com_arr_y.push_back(arr_even_y[i] - arr_odd_y[i]);
    }
    /*for (int i = 0; i < com_arr_y.size(); i++)
    {
        cout << arr_even_x[i] << "," << com_arr_y[i] << endl;
    }*/
}

void imag_part_com_ip_qp(vector<float> &arr_even_x, vector<float> &arr_even_y_imag, vector<float> &arr_odd_y_imag) {

    /*** Values of the Q(t)Cos(wct) and I(t)Sin(wct) computed from the prior functions are added to
    compute the combined imaginary part of the passband modulated waveform in time domain. ***/
    for(int i=0; i<arr_even_x.size(); i++) {
        com_arr_y_imag.push_back(arr_even_y_imag[i] + arr_odd_y_imag[i]);
    }

    /*for (int i=0;i < arr_even_x.size();i++){
        cout << arr_even_x[i]  << ',' << com_arr_y_imag[i] << endl;
    }*/
}

/*** GNUplot Function to plot the generated waveforms. (Issues with generating plots for large number of bits/symbols)
 Shifted to the MATLAB plotting from csv file having the values generated from this C++ code. ***/

void plotResults(vector<float> &xData, vector<float> &yData, char *tempDataFileName)
{
    FILE *gnuplotPipe, *tempDataFile;
    double x, y;
    int i;
    gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe)
    {
        fprintf(gnuplotPipe, "plot \"%s\" with lines\n", tempDataFileName);
        fflush(gnuplotPipe);
        tempDataFile = fopen(tempDataFileName, "w");
        for (i = 0; i < xData.size(); i++)
        {
            x = xData[i];
            y = yData[i];
            fprintf(tempDataFile, "%lf %lf\n", x, y);
            cout << x << ',' << y << endl;
        }
        fclose(tempDataFile);
        printf("press enter to continue...");
        getchar();
        remove(tempDataFileName);
        fprintf(gnuplotPipe, "exit \n"); 
    }
    else
    {
        printf("gnuplot not found...");
    }
}

int main()
{
    /*** Generation of random input bit sequence to simulate PPDU bits ***/
    
    srand(1234);  // Set 1234 in the random seed for reproducibility of results.
    int randn = 0;
    int num_bits = 2048;
    int randns[2048] = {};
    
    for (int i = 0; i < num_bits; i++)
    {
        randn = (rand() >> i) & 1;
        randns[i] = randn;
    }

    /*** Declaration of the vectors for storing the x,y values of the baseband, passband inphase and quadrature components 
     for multiple symbols. ***/

    vector<float> mul_arr_even_x;
    vector<float> mul_arr_even_y; 
    vector<float> mul_arr_odd_x; 
    vector<float> mul_arr_odd_y; 
    vector<float> mul_arr_even_y_imag;
    vector<float> mul_arr_odd_y_imag;

    ofstream fout;
    fout.open("test_spc.csv");
    
    int syms = (sizeof(randns)/sizeof(int))/4;  //Number of symbols from bits.

    /*** Setting the values of the parameters Tc (Chip period), fc (Carrier Frequency), spc(samplesPerChip) ***/

    float Tc = 0.5;   // in microseconds.
    long long fc = 2.45 * 1e9;  // in Hz.
    int spc = 8;

    /*** Defining the variables for tracking the count variable for correct time t according to the number of symbols. ***/
    
    int curr_count_even = 0;
    int curr_count_odd = 0;

    // Generation of OQPSK baseband and passband plots for multiple symbols. 

    for(int i=0;i<syms;i++){
        vector<int> data_sym;
        for (int j=0;j<4;j++){   // Each data symbol contains 4 bits.
            data_sym.push_back(randns[j+4*i]);
        }
        int count_even = curr_count_even;
        int count_odd = curr_count_odd;

        /*** Call to the bit to symbol map function to get symbol. ***/

        int val = btsm(data_sym);

        /*** Call to the symbol to chip map function to get the corresponding chip sequences. ***/

        string chip_map = stcm(val);

        /*** Call to the function to split the chip sequence into even and odd numbered chip sequences. ***/

        Struct e_o_string;
        e_o_string = c_seq_split(chip_map);  

        /*** Call to the functions for OQPSK modulated baseband waveform generation. Comment when running passband case. ***/
        curr_count_even = oqpsk_ip(e_o_string.even,Tc, count_even, spc);
        curr_count_odd = oqpsk_qp(e_o_string.odd,Tc, count_odd, spc);

        /***Inserting the values in the vector for multiple symbols by concatenating the vectors from the end of 
         the array for multiple symbols array and starting of the correpsonding single symbol vector. ***/

        mul_arr_even_x.insert(mul_arr_even_x.end(), arr_even_x.begin(), arr_even_x.end());
        mul_arr_even_y.insert(mul_arr_even_y.end(), arr_even_y.begin(), arr_even_y.end());
        mul_arr_odd_x.insert(mul_arr_odd_x.end(), arr_odd_x.begin(), arr_odd_x.end());
        mul_arr_odd_y.insert(mul_arr_odd_y.end(), arr_odd_y.begin(), arr_odd_y.end());

        /*** Call to the functions for OQPSK modulated passband waveform generation. Comment when running baseband case. ***/

        /*curr_count_even = com_ip_oqpsk(e_o_string.even, Tc, fc, count_even, spc);
        float step_size = Tc/spc;
        if(i == syms-1){
            // Padding zeroes at the end of modulated inphase component of oqpsk for multiple symbols
            // to facilitate the subtraction of the inphase and quadrature phase components to get the 
            // real part of the passband modulated RF signal.
            for (int i = 2 * curr_count_even * spc + 1; i <= 2 * curr_count_even * spc + spc; i++) {
                arr_even_x.push_back(i * step_size);
                arr_even_y.push_back(0);
                arr_even_y_imag.push_back(0);
            }
        }
        mul_arr_even_x.insert(mul_arr_even_x.end(), arr_even_x.begin(), arr_even_x.end());
        mul_arr_even_y.insert(mul_arr_even_y.end(), arr_even_y.begin(), arr_even_y.end());
        mul_arr_even_y_imag.insert(mul_arr_even_y_imag.end(), arr_even_y_imag.begin(), arr_even_y_imag.end());
        curr_count_odd = com_qp_oqpsk(e_o_string.odd, Tc, fc, count_odd, spc);
        mul_arr_odd_x.insert(mul_arr_odd_x.end(), arr_odd_x.begin(), arr_odd_x.end());
        mul_arr_odd_y.insert(mul_arr_odd_y.end(), arr_odd_y.begin(), arr_odd_y.end());
        mul_arr_odd_y_imag.insert(mul_arr_odd_y_imag.end(), arr_odd_y_imag.begin(), arr_odd_y_imag.end());*/

        /*** Clearing the single symbol vectors to avoid duplication while inserting the values for multiple symbols. ***/

        arr_even_x.clear();
        arr_even_y.clear();
        arr_odd_x.clear();
        arr_odd_y.clear();
        arr_even_y_imag.clear();
        arr_odd_y_imag.clear();
    }
    
    /*** Call to function to combine the values of modulated I(t) and Q(t) to form the real part of passband signal. ***/
    //real_part_com_ip_qp(mul_arr_even_x, mul_arr_even_y, mul_arr_odd_y);
    
    /*** Call to function to combine the values of modulated I(t) and Q(t) to form the imaginary part of passband signal. ***/
    //imag_part_com_ip_qp(mul_arr_even_x, mul_arr_even_y_imag, mul_arr_odd_y_imag);
    
    /*** Iterative loop for the exporting of values onto a csv file for plotting in MATLAB ***/

    for (int i=0;i<mul_arr_odd_x.size();i++){
        fout << mul_arr_odd_x[i] << "," << mul_arr_even_y[i] << "," << mul_arr_odd_y[i] << endl;
    }

    /*** Call to the GNU plot functions for plotting the generated OQPSK waveforms. (For Initial Plotting) ***/

    //plotResults(mul_arr_even_x, mul_arr_even_y, "Inphase Component");
    //plotResults(mul_arr_odd_x, mul_arr_odd_y, "Quadrature Component");
    //plotResults(mul_arr_even_x, com_arr_y, "Real part of combined modulated signal");

    for(int i=0;i<mul_arr_odd_x.size();i++){
        cout << mul_arr_odd_x[i] <<  "," << mul_arr_odd_y[i] << "," << mul_arr_even_y[i] << endl; 
    }
    
    //cout << mul_arr_odd_x.size() << endl;

    return 0;
}