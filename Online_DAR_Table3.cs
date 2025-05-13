using System;
using System.Collections.Generic;
using System.Xml;
using System.Text;
using System.IO;
using System.Data;
using System.Data.OleDb;
using ILOG.Concert;
using ILOG.CPLEX;
using System.Runtime.InteropServices;
using System.Security.Cryptography.X509Certificates;
using System.Security.Cryptography;
using System.Data.SqlTypes;
using System.Runtime.Serialization.Formatters;
//=======================
using System.Threading.Tasks;

using Itinero;
using Itinero.IO.Osm;
using Itinero.Osm.Vehicles;
using System.Dynamic;
using Itinero.LocalGeo;
using System.Threading;
using Itinero.Algorithms.Networks.Analytics.Heatmaps;
using System.Management.Instrumentation;

namespace Cplex_Integrated_BAP_SAP
{

    class Program
    {
        public const int W = 1000; 
        public const short n = 2000; 

        public const short nprime = 10;  
        public const short Q = 4;  
        public const double apha = 0.5; 
        public const int discre_no = 10;

        public const double repo_or_not_par = 0.005; 
        public const double T_0 = 1.5 * 3600;

        public int[,] selection_j = new int[n + 1, 50]; 

        public int[] true_j = new int[ n + 1]; 

        public double[] aao_j = new double[ n + 1];
        public double[] bbo_j = new double[ n + 1];
        public double[] aad_j = new double[ n + 1];
        public double[] bbd_j = new double[ n + 1];

        public int[] q_j = new int[2 * (n + 1)]; 

        public double[] sw_j = new double[n + 1]; 
        public double[] lw_j = new double[2 * (n + 1)]; 
        public double[] T_j = new double[n + 1]; 
        public double[] g_j = new double[n + 1]; 
        public double[,] t_ij = new double[2 * (n + 1), 2 * (n + 1)];  
        public double[,] pi_ij = new double[2 * (n + 1), 2 * (n + 1)]; 
        public double[,] pipi_ij = new double[2 * (n + 1), 2 * (n + 1)]; 

        public const double cost = 3.2 / 3600.0;


        public double lwj_value = 10000000;

        public const int L = (int)268435400 / nwords / 2;
        public int[] prel_l = new int[L];
        public short[] d_l = new short[L];


        public const int nwords = (int)(n / 30) + 1;
        public int nwords_1 = (int)(Math.Ceiling(W / 30.0) - 1);
        public int nwords_2 = 0;


        public int[,] preset1_ll = new int[L, nwords + 1];
        public int[,] preset2_ll = new int[L, nwords + 1];
        public int[,] preset3_ll = new int[L, nwords + 1];

        public int[] n_l = new int[L];
        public int[] Lall = new int[1];
        public double[] T_l = new double[L];
        public double[] pipipi_l = new double[L];
        public int[] S_l = new int[L];
        public int[] Q_l = new int[L];

        public double[] aa_ln = new double[100];
        public double[] bb_ln = new double[100];
        public double[] case_ln = new double[100];

        public int previous_drop_off_location;
        public int previous_drop_off_point;
        public double T_arriving_previous_drop_off;
        public double g_arriving_previous_drop_off;

        public int[] lbit_j = new int[n + 1];
        public int[] lwrd_j = new int[n + 1];


        public string strConnection;


        public double[,] candidate_aam_jc = new double[2 * n + 2, discre_no + 1];
        public double[,] candidate_bbm_jc = new double[2 * n + 2, discre_no + 1];

        public TextWriter writer = File.CreateText("output.txt"); 

        public int ride_sharing = -1;
        public int Fle_pickup = -1;
        public int Fle_dropoff = -1;
        public int repositioning = -1;
        public void overall_control(int all_control)
        {
            if (all_control == 1)
            {
                ride_sharing = 0;
                Fle_pickup = 0;
                Fle_dropoff = 0;
                repositioning = 1;
            }

            if (all_control == 2)
            {
                ride_sharing = 0;
                Fle_pickup = 1;
                Fle_dropoff = 1;
                repositioning = 1;
            }

            if (all_control == 3)
            {
                ride_sharing = 1;
                Fle_pickup = 0;
                Fle_dropoff = 0;
                repositioning = 1;
            }

            if (all_control == 4)
            {
                ride_sharing = 1;
                Fle_pickup = 1;
                Fle_dropoff = 1;
                repositioning = 1;
            }

            if (all_control == 5)
            {
                ride_sharing = 1;
                Fle_pickup = 1;
                Fle_dropoff = 0;
                repositioning = 1;
            }

            if (all_control == 6)
            {
                ride_sharing = 1;
                Fle_pickup = 0;
                Fle_dropoff = 1;
                repositioning = 1;
            }

            if (all_control == 7)
            {
                ride_sharing = 1;
                Fle_pickup = 0;
                Fle_dropoff = 0;
                repositioning = 0;
            }

            if (all_control == 8)
            {
                ride_sharing = 1;
                Fle_pickup = 1;
                Fle_dropoff = 1;
                repositioning = 0;
            }

            if (all_control == 9)
            {
                ride_sharing = 0;
                Fle_pickup = 1;
                Fle_dropoff = 0;
                repositioning = 1;
            }

            if (all_control == 10)
            {
                ride_sharing = 0;
                Fle_pickup = 0;
                Fle_dropoff = 1;
                repositioning = 1;
            }

            if (all_control == 11)
            {
                ride_sharing = 0;
                Fle_pickup = 0;
                Fle_dropoff = 0;
                repositioning = 0;
            }

            if (all_control == 12)
            {
                ride_sharing = 0;
                Fle_pickup = 1;
                Fle_dropoff = 1;
                repositioning = 0;
            }

            if (all_control == 13)
            {
                ride_sharing = 1;
                Fle_pickup = 1;
                Fle_dropoff = 0;
                repositioning = 0;
            }

            if (all_control == 14)
            {
                ride_sharing = 1;
                Fle_pickup = 0;
                Fle_dropoff = 1;
                repositioning = 0;
            }

            if (all_control == 15)
            {
                ride_sharing = 0;
                Fle_pickup = 1;
                Fle_dropoff = 0;
                repositioning = 0;
            }

            if (all_control == 16)
            {
                ride_sharing = 0;
                Fle_pickup = 0;
                Fle_dropoff = 1;
                repositioning = 0;
            }

            Console.WriteLine("all_control = " + all_control);
            Console.WriteLine("ride_sharing = " + ride_sharing + "  Fle_pickup = " + Fle_pickup + "  Fle_dropoff = " + Fle_dropoff + "  repositioning = " + repositioning);
        }
    
        double Real_Coordinate_Time_Car_data(int jj_1, int jj_1_candidate_no, int jj_2, int jj_2_candidate_no)
        {
            int intersection_1;
            int intersection_2;

            if (jj_1 <= n)
            {
                if (true_j[jj_1] > 100000)
                    intersection_1 = true_j[jj_1] - 100000;
                else
                    intersection_1 = pickupi_match_j_no[true_j[jj_1], jj_1_candidate_no];

            }
            else
            {
                if (true_j[jj_1 - n] > 100000)
                    intersection_1 = true_j[jj_1 - n] - 100000;
                else
                    intersection_1 = dropoffi_match_j_no[true_j[jj_1 - n], jj_1_candidate_no];
            }


            if (jj_2 <= n)
            {
                if (true_j[jj_2] > 100000)
                    intersection_2 = true_j[jj_2] - 100000;
                else
                    intersection_2 = pickupi_match_j_no[true_j[jj_2], jj_2_candidate_no];
            }
            else
            {
                if (true_j[jj_2 - n] > 100000)
                    intersection_2 = true_j[jj_2 - n] - 100000;
                else
                    intersection_2 = dropoffi_match_j_no[true_j[jj_2 - n], jj_2_candidate_no];
            }

            return Intersection_car_t_ij[intersection_1, intersection_2];
        }

        double Real_Coordinate_Dis_Car_data(int jj_1, int jj_1_candidate_no, int jj_2, int jj_2_candidate_no)
        {
            int intersection_1;
            int intersection_2;
            if (jj_1 <= n)
            {
                if (true_j[jj_1] > 100000)
                    intersection_1 = true_j[jj_1] - 100000;
                else
                    intersection_1 = pickupi_match_j_no[true_j[jj_1], jj_1_candidate_no];

            }
            else
            {
                if (true_j[jj_1 - n] > 100000)
                    intersection_1 = true_j[jj_1 - n] - 100000;
                else
                    intersection_1 = dropoffi_match_j_no[true_j[jj_1 - n], jj_1_candidate_no];
            }


            if (jj_2 <= n)
            {
                if (true_j[jj_2] > 100000)
                    intersection_2 = true_j[jj_2] - 100000;
                else
                    intersection_2 = pickupi_match_j_no[true_j[jj_2], jj_2_candidate_no];
            }
            else
            {
                if (true_j[jj_2 - n] > 100000)
                    intersection_2 = true_j[jj_2 - n] - 100000;
                else
                    intersection_2 = dropoffi_match_j_no[true_j[jj_2 - n], jj_2_candidate_no];
            }

            return Intersection_car_d_ij[intersection_1, intersection_2];
        }

        double Real_Coordinate_Dis_Car_data_j_truej(int jj_1, int jj_1_candidate_no, int truejj_2, int truejj_2_candidate_no)
        {
            int intersection_1;
            int intersection_2;

            if (jj_1 <= n)
            {
                if (true_j[jj_1] > 100000)
                    intersection_1 = true_j[jj_1] - 100000;
                else
                    intersection_1 = pickupi_match_j_no[true_j[jj_1], jj_1_candidate_no];

            }
            else
            {
                if (true_j[jj_1 - n] > 100000)
                    intersection_1 = true_j[jj_1 - n] - 100000;
                else
                    intersection_1 = dropoffi_match_j_no[true_j[jj_1 - n], jj_1_candidate_no];

            }

            if (truejj_2 <= truen)
            {
                intersection_2 = pickupi_match_j_no[truejj_2, truejj_2_candidate_no];
            }
            else if (truejj_2 > truen && truejj_2 <= 100000)
            {
                intersection_2 = dropoffi_match_j_no[truejj_2 - truen, truejj_2_candidate_no];
            }
            else
            {
                intersection_2 = truejj_2 - 100000;
            }

            return Intersection_car_d_ij[intersection_1, intersection_2];
        }

        double Real_Coordinate_Dis_Car_data_intersection_truej(int intersection, int truejj_2, int truejj_2_candidate_no)
        {
            int intersection_1;
            int intersection_2;

            intersection_1 = intersection;

            if (truejj_2 <= truen)
            {
                intersection_2 = pickupi_match_j_no[truejj_2, truejj_2_candidate_no];
            }
            else if (truejj_2 > truen && truejj_2 <= 100000)
            {
                intersection_2 = dropoffi_match_j_no[truejj_2 - truen, truejj_2_candidate_no];
            }
            else
            {
                intersection_2 = truejj_2 - 100000;
            }

            return Intersection_car_d_ij[intersection_1, intersection_2];
        }

        double Real_Coordinate_Time_Car_data_j_intersection(int jj_1, int jj_1_candidate_no, int intersection)
        {
            int intersection_1;
            int intersection_2;

            if (jj_1 <= n)
            {
                if (true_j[jj_1] > 100000)
                    intersection_1 = true_j[jj_1] - 100000;
                else
                    intersection_1 = pickupi_match_j_no[true_j[jj_1], jj_1_candidate_no];
            }
            else
            {
                if (true_j[jj_1 - n] > 100000)
                    intersection_1 = true_j[jj_1 - n] - 100000;
                else
                    intersection_1 = dropoffi_match_j_no[true_j[jj_1 - n], jj_1_candidate_no];
            }

            intersection_2 = intersection;

            return Intersection_car_t_ij[intersection_1, intersection_2];
        }

        double Real_Coordinate_Time_Pedestrian_data(int jj, int jj_candidate_no)
        {
            int intersection_1;
            int intersection_2;

            if (jj <= n)
            {
                intersection_1 = pickupi_match_j_no[true_j[jj], 1];
            }
            else
            {
                intersection_1 = dropoffi_match_j_no[true_j[jj - n], 1];
            }


            if (jj <= n)
            {
                intersection_2 = pickupi_match_j_no[true_j[jj], jj_candidate_no];
            }
            else
            {
                intersection_2 = dropoffi_match_j_no[true_j[jj - n], jj_candidate_no];
            }

            if (jj <= n)
            {
                return Intersection_ped_t_ij[intersection_1, intersection_2];
            }
            else
            {
                return Intersection_ped_t_ij[intersection_2, intersection_1];
            }
        }

        double Real_Coordinate_Dis_Pedestrian_data(int jj, int jj_candidate_no)
        {
            int intersection_1;
            int intersection_2;

            if (jj <= n)
            {
                intersection_1 = pickupi_match_j_no[true_j[jj], 1];
            }
            else
            {
                intersection_1 = dropoffi_match_j_no[true_j[jj - n], 1];
            }


            if (jj <= n)
            {
                intersection_2 = pickupi_match_j_no[true_j[jj], jj_candidate_no];
            }
            else
            {
                intersection_2 = dropoffi_match_j_no[true_j[jj - n], jj_candidate_no];
            }

            if (jj <= n)
            {
                return Intersection_ped_d_ij[intersection_1, intersection_2];
            }
            else
            {
                return Intersection_ped_d_ij[intersection_2, intersection_1];
            }
        }

        double Real_Coordinate_Dis_Pedestrian_data_truej(int truej, int truej_candidate_no)
        {
            int intersection_1;
            int intersection_2;

            if (truej <= truen)
            {
                intersection_1 = pickupi_match_j_no[truej, 1];
            }
            else
            {
                intersection_1 = dropoffi_match_j_no[truej - truen, 1];
            }


            if (truej <= truen)
            {
                intersection_2 = pickupi_match_j_no[truej, truej_candidate_no];
            }
            else
            {
                intersection_2 = dropoffi_match_j_no[truej - truen, truej_candidate_no];
            }

            if (truej <= truen)
            {
                return Intersection_ped_d_ij[intersection_1, intersection_2];
            }
            else
            {
                return Intersection_ped_d_ij[intersection_2, intersection_1];
            }
        }


        public int check_opt = 0;
  
 
        public const int R = 67000; 
        public double[] Requesttime_j = new double[n + 1]; 
        public double[] Deadline_j = new double[n + 1]; 
        public int[] earliest_time_j = new int[n + 1];
        public int[] latest_time_j = new int[n + 1]; 
        public int r = 0;
        public int[] start_r = new int[R + 1]; 
        public int[] end_r = new int[R + 1]; 
        public double[] g_r = new double[R + 1]; 
        public int r_dummy = 0;
        public int[] start_r_dummy = new int[R + 1]; 
        public int[] end_r_dummy = new int[R + 1]; 
        public double[] g_r_dummy = new double[R + 1]; 
        public int[,] a_rj = new int[R + 1, 2 * n + 2]; 
        public int[,] seq_rj = new int[R + 1, 2 * n + 2]; 
        public int[] record_r = new int[R + 1]; 
        public float[,] T_rj = new float[R + 1, 2 * n + 2]; 
        public float[,] pipipi_rj = new float[R + 1, 2 * n + 2]; 


        public const int S = 100000;
        public int s = 0;
        public int[] n_location_s = new int[S + 1];
        public int[] time_s = new int[S + 1];
        public int[] check_s = new int[S + 1];
        public const double unit = 10;


        public int[] set_s = new int[S + 1];
        public int[] vehicle_s = new int[W + 1];

 
        public int d1 = 0;
        public int last_d1 = 0;

    
        public void input_vehicle_origin_at_intersections(int no_vehicle)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            strConnection += @"Data Source=|DataDirectory|\Manhattan_Intersection_Coordinates_Mod.mdb";

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            d1 = 0;
            Random rnd = new Random(123);
            int rnd1;
            for (int vv = 1; vv < no_vehicle + 1; vv++)
            {
                d1 = d1 + 1;

                double rnd_seed = rnd.NextDouble();
                rnd1 = (int)Math.Ceiling((0) * rnd_seed + (Total_I) * (1 - rnd_seed));

                string str1 = "select * from Manhattan_Intersection_Coordinates";
                OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
                OleDbDataReader myReader001 = myCommand001.ExecuteReader();

                while (myReader001.Read())
                {
                    if (myReader001.GetInt16(1) == rnd1)
                    {
                        status_j[d1] = 100;
                        true_j[d1] = 100000 + rnd1;

                        truej_vehicle_stop[vv, 0] = true_j[d1];

                        aao_j[d1] = myReader001.GetDouble(2); 
                        bbo_j[d1] = myReader001.GetDouble(3);
                        aad_j[d1] = myReader001.GetDouble(2); 
                        bbd_j[d1] = myReader001.GetDouble(3); 

                        break;
                    }
                }
                myReader001.Close();
            }
            objConnection.Close();
        }

        public const int I = 4000;
        public int Total_I = 0;
        public double[] aac_i = new double[I];
        public double[] bbc_i = new double[I];
        public void input_intersections()
        {

            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            strConnection += @"Data Source=|DataDirectory|\Manhattan_Intersection_Coordinates_Mod.mdb";

            
            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Manhattan_Intersection_Coordinates"; 
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();

            int dd1 = 0;
            while (myReader001.Read())
            {
                dd1 = myReader001.GetInt16(1);
                aac_i[dd1] = myReader001.GetDouble(2);
                bbc_i[dd1] = myReader001.GetDouble(3);
                if (Total_I < dd1)
                    Total_I = dd1; 

            }
            myReader001.Close();
            objConnection.Close();
        }

        public double[,] Intersection_car_t_ij = new double[I, I]; 
        public double[,] Intersection_car_d_ij = new double[I, I];
        public double[,] Intersection_ped_t_ij = new double[I, I];
        public double[,] Intersection_ped_d_ij = new double[I, I];

        public int[,] pickupi_match_j_no = new int[truen + 1, discre_no + 1];
        public int[,] dropoffi_match_j_no = new int[truen + 1, discre_no + 1];

        public void input_ball_Discre0611_justinput_rolling_coordinates_matching_changej(double previous_resolve_start, double resolve_start, int no_vehicle, int Fle_pickup, int Fle_dropoff)
        {
            last_d1 = d1;
            d1 = 0;
            for (int j1 = 1; j1 < last_d1 + 1; j1++)
            {
                 if (status_j[j1] == 100)
                 {
                    d1 = d1 + 1;

                 
                    for (int v = 1; v < no_vehicle + 1; v++)
                    {
                        if (n_location_v[v] == j1 + n)
                        {
                            n_location_v[v] = d1 + n;
                        }
                    }
           

                    status_j[d1] = status_j[j1];
                    true_j[d1] = true_j[j1];
                    Requesttime_j[d1] = Requesttime_j[j1];
                    q_j[d1] = q_j[j1];
                    g_j[d1] = g_j[j1];

                    aao_j[d1] = aao_j[j1];
                    bbo_j[d1] = bbo_j[j1];
                    aad_j[d1] = aad_j[j1];
                    bbd_j[d1] = bbd_j[j1];

                    for (int cc = 1; cc < 10; cc++)
                    {
                        candidate_aam_jc[d1, cc] = candidate_aam_jc[j1, cc];
                        candidate_bbm_jc[d1, cc] = candidate_bbm_jc[j1, cc];
                        candidate_aam_jc[d1 + n, cc] = candidate_aam_jc[j1 + n, cc];
                        candidate_bbm_jc[d1 + n, cc] = candidate_bbm_jc[j1 + n, cc];

                    }
                    candidate_aam_jc[d1, 10] = aao_j[j1];
                    candidate_bbm_jc[d1, 10] = bbo_j[j1];
                    candidate_aam_jc[d1 + n, 10] = aad_j[j1];
                    candidate_bbm_jc[d1 + n, 10] = bbd_j[j1];

                    sw_j[d1] = 4.0 / 3600.0;

   
                    if (reposition_j[j1] == -1)
                    {
                        Requesttime_j[d1] = 0;
                        q_j[d1] = 0;
                        g_j[d1] = 0;
                        for (int cc = 1; cc < 10; cc++)
                        {
                            candidate_aam_jc[d1, cc] = 0;
                            candidate_bbm_jc[d1, cc] = 0;
                            candidate_aam_jc[d1 + n, cc] = 0;
                            candidate_bbm_jc[d1 + n, cc] = 0;
                        }
                    }
                    reposition_j[d1] = 0;
                }
            }

            int[] update_v = new int[W + 1];
            for (int vv = 1; vv < W + 1; vv++)
            {
                update_v[vv] = -1;

            }
            for (int j1 = 1; j1 < last_d1 + 1; j1++)
            {
                if (status_j[j1] == 88)
                {
                    d1 = d1 + 1;

                    status_j[d1] = status_j[j1];

                    cust_set_j[d1] = cust_set_j[j1];

                    for (int vv = 1; vv < W + 1; vv++)
                    {
                        if (first_cust_v[vv] == j1 && update_v[vv] == -1 && first_cust_v[vv] != 0)
                        {
                            first_cust_v[vv] = d1;
                            update_v[vv] = 1;
                            break;
                        }
                    }

                    true_j[d1] = true_j[j1];
                    Requesttime_j[d1] = Requesttime_j[j1];
                    q_j[d1] = q_j[j1];
                    g_j[d1] = g_j[j1];

                    aao_j[d1] = aao_j[j1];
                    bbo_j[d1] = bbo_j[j1];
                    aad_j[d1] = aad_j[j1];
                    bbd_j[d1] = bbd_j[j1];

                    for (int cc = 1; cc < 10; cc++)
                    {
                        candidate_aam_jc[d1, cc] = candidate_aam_jc[j1, cc];
                        candidate_bbm_jc[d1, cc] = candidate_bbm_jc[j1, cc];
                        candidate_aam_jc[d1 + n, cc] = candidate_aam_jc[j1 + n, cc];
                        candidate_bbm_jc[d1 + n, cc] = candidate_bbm_jc[j1 + n, cc];

                    }
                    candidate_aam_jc[d1, 10] = aao_j[j1];
                    candidate_bbm_jc[d1, 10] = bbo_j[j1];
                    candidate_aam_jc[d1 + n, 10] = aad_j[j1];
                    candidate_bbm_jc[d1 + n, 10] = bbd_j[j1];

                    sw_j[d1] = 4.0 / 3600.0;

                   
                }
            }

          
            for (int jjj = d1 + 1; jjj < n + 1; jjj++)
            {
                status_j[jjj] = 0;
            }
        
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";

        
            strConnection += @"Data Source=|DataDirectory|\Instance_file_Mod.mdb";


            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Instance_file"; 
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();

            while (myReader001.Read())
            {
                if (myReader001.GetDouble(2) <= resolve_start && myReader001.GetDouble(2) > previous_resolve_start)
                {
                    d1 = d1 + 1; 

                    status_j[d1] = 55;

                    cust_set_j[d1] = 0;

                    true_j[d1] = myReader001.GetInt16(1);

                    Requesttime_j[d1] = myReader001.GetDouble(2);
                    q_j[d1] = myReader001.GetInt16(3);
                    g_j[d1] = myReader001.GetDouble(4);

                    aao_j[d1] = myReader001.GetDouble(5);
                    bbo_j[d1] = myReader001.GetDouble(6);
                    aad_j[d1] = myReader001.GetDouble(7);
                    bbd_j[d1] = myReader001.GetDouble(8);

     
                    max_truej_untilnow = true_j[d1];
                    aao_truej[true_j[d1]] = aao_j[d1];
                    bbo_truej[true_j[d1]] = bbo_j[d1];
                    aad_truej[true_j[d1]] = aad_j[d1];
                    bbd_truej[true_j[d1]] = bbd_j[d1];
      

                    candidate_aam_jc[d1, 1] = myReader001.GetDouble(9);
                    candidate_bbm_jc[d1, 1] = myReader001.GetDouble(10);

                    candidate_aam_jc[d1, 2] = myReader001.GetDouble(11);
                    candidate_bbm_jc[d1, 2] = myReader001.GetDouble(12);

                    candidate_aam_jc[d1, 3] = myReader001.GetDouble(13);
                    candidate_bbm_jc[d1, 3] = myReader001.GetDouble(14);

                    candidate_aam_jc[d1, 4] = myReader001.GetDouble(15);
                    candidate_bbm_jc[d1, 4] = myReader001.GetDouble(16);

                    candidate_aam_jc[d1, 5] = myReader001.GetDouble(17);
                    candidate_bbm_jc[d1, 5] = myReader001.GetDouble(18);

                    candidate_aam_jc[d1, 6] = myReader001.GetDouble(19);
                    candidate_bbm_jc[d1, 6] = myReader001.GetDouble(20);

                    candidate_aam_jc[d1, 7] = myReader001.GetDouble(21);
                    candidate_bbm_jc[d1, 7] = myReader001.GetDouble(22);

                    candidate_aam_jc[d1, 8] = myReader001.GetDouble(23);
                    candidate_bbm_jc[d1, 8] = myReader001.GetDouble(24);

                    candidate_aam_jc[d1, 9] = myReader001.GetDouble(25);
                    candidate_bbm_jc[d1, 9] = myReader001.GetDouble(26);

                    candidate_aam_jc[d1 + n, 1] = myReader001.GetDouble(27);
                    candidate_bbm_jc[d1 + n, 1] = myReader001.GetDouble(28);

                    candidate_aam_jc[d1 + n, 2] = myReader001.GetDouble(29);
                    candidate_bbm_jc[d1 + n, 2] = myReader001.GetDouble(30);

                    candidate_aam_jc[d1 + n, 3] = myReader001.GetDouble(31);
                    candidate_bbm_jc[d1 + n, 3] = myReader001.GetDouble(32);

                    candidate_aam_jc[d1 + n, 4] = myReader001.GetDouble(33);
                    candidate_bbm_jc[d1 + n, 4] = myReader001.GetDouble(34);

                    candidate_aam_jc[d1 + n, 5] = myReader001.GetDouble(35);
                    candidate_bbm_jc[d1 + n, 5] = myReader001.GetDouble(36);

                    candidate_aam_jc[d1 + n, 6] = myReader001.GetDouble(37);
                    candidate_bbm_jc[d1 + n, 6] = myReader001.GetDouble(38);

                    candidate_aam_jc[d1 + n, 7] = myReader001.GetDouble(39);
                    candidate_bbm_jc[d1 + n, 7] = myReader001.GetDouble(40);

                    candidate_aam_jc[d1 + n, 8] = myReader001.GetDouble(41);
                    candidate_bbm_jc[d1 + n, 8] = myReader001.GetDouble(42);

                    candidate_aam_jc[d1 + n, 9] = myReader001.GetDouble(43);
                    candidate_bbm_jc[d1 + n, 9] = myReader001.GetDouble(44);

                    candidate_aam_jc[d1, 10] = aao_j[d1];
                    candidate_bbm_jc[d1, 10] = bbo_j[d1];
                    candidate_aam_jc[d1 + n, 10] = aad_j[d1];
                    candidate_bbm_jc[d1 + n, 10] = bbd_j[d1];

                    sw_j[d1] = 4.0 / 3600.0; 
 
                }

                if (myReader001.GetDouble(2) > resolve_start)
                    break;
            }

            myReader001.Close();
            objConnection.Close();

            aao_j[0] = -73.9464428326014;
            bbo_j[0] = 40.7818193234921;
            aad_j[0] = -73.9464428326014;
            bbd_j[0] = 40.7818193234921;

            lwj_value = 1000000;

            if (Fle_pickup == 1 && Fle_dropoff == 1)
            {
                for (int j = 1; j < d1 + 1; j++)
                {
                    lw_j[j] = lwj_value;
                }

                for (int j = n + 1; j < 2 * n + 1; j++)
                {
                    lw_j[j] = lwj_value;
                }
            }
            if (Fle_pickup == 1 && Fle_dropoff == 0)
            {
                for (int j = 1; j < d1 + 1; j++)
                {
                    lw_j[j] = lwj_value;
                }

                for (int j = n + 1; j < 2 * n + 1; j++)
                {
                    lw_j[j] = 0;
                }
            }
            if (Fle_pickup == 0 && Fle_dropoff == 1)
            {
                for (int j = 1; j < d1 + 1; j++)
                {
                    lw_j[j] = 0;
                }

                for (int j = n + 1; j < 2 * n + 1; j++)
                {
                    lw_j[j] = lwj_value;
                }
            }
            if (Fle_pickup == 0 && Fle_dropoff == 0)
            {
                for (int j = 1; j < d1 + 1; j++)
                {
                    lw_j[j] = 0;
                }

                for (int j = n + 1; j < 2 * n + 1; j++)
                {
                    lw_j[j] = 0;
                }
            }

            for (int j = 1; j < d1 + 1; j++)
            {
                T_j[j] = Real_Coordinate_Time_Car_data(j, 1, j + n, 1);
            }

            for (int i = 1; i < d1 + 1; i++)
            {
                if (i < n + 1)
                {
                    t_ij[i, 2 * n + 1] = 0;

                    pi_ij[i, 2 * n + 1] = t_ij[i, 2 * n + 1] - g_j[i] / 2;

                    t_ij[0, i] = 0;

                    pi_ij[0, i] = t_ij[0, i] - g_j[i] / 2;
                }

                if (i >= n + 1)
                {
                    t_ij[i, 2 * n + 1] = 0;

                    pi_ij[i, 2 * n + 1] = t_ij[i, 2 * n + 1];

                    t_ij[0, i] = 0;

                    pi_ij[0, i] = t_ij[0, i];
                }

                for (int j = 1; j < d1 + 1; j++)
                {
                    if (i < n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2 - g_j[j] / 2;
                    }

                    if (i >= n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[j] / 2;
                    }
                    if (i < n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2;
                    }
                    if (i >= n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j];
                    }
                }

                for (int j = n + 1; j < (n + d1) + 1; j++)
                {
                    if (i < n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2 - g_j[j] / 2;
                    }

                    if (i >= n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[j] / 2;
                    }
                    if (i < n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2;
                    }
                    if (i >= n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j];
                    }
 
                }

            }

            for (int i = n + 1; i < (n + d1) + 1; i++)
            {
                if (i < n + 1)
                {
                    t_ij[i, 2 * n + 1] = 0;

                    pi_ij[i, 2 * n + 1] = t_ij[i, 2 * n + 1] - g_j[i] / 2;

                    t_ij[0, i] = 0;

                    pi_ij[0, i] = t_ij[0, i] - g_j[i] / 2;
                }

                if (i >= n + 1)
                {
                    t_ij[i, 2 * n + 1] = 0;

                    pi_ij[i, 2 * n + 1] = t_ij[i, 2 * n + 1];

                    t_ij[0, i] = 0;

                    pi_ij[0, i] = t_ij[0, i];
                }

                for (int j = 1; j < d1 + 1; j++)
                {
                    if (i < n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2 - g_j[j] / 2;
                    }

                    if (i >= n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[j] / 2;
                    }
                    if (i < n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2;
                    }
                    if (i >= n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j];
                    }
                }

                for (int j = n + 1; j < (n + d1) + 1; j++)
                {
                    if (i < n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2 - g_j[j] / 2;
                    }

                    if (i >= n + 1 && j < n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[j] / 2;
                    }
                    if (i < n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j] - g_j[i] / 2;
                    }
                    if (i >= n + 1 && j >= n + 1)
                    {
                        t_ij[i, j] = Real_Coordinate_Time_Car_data(i, 1, j, 1);

                        pi_ij[i, j] = t_ij[i, j];
                    }

                }
            }

            for (int j = n + 1; j < 2 * n + 1; j++)
            {
                q_j[j] = -q_j[j - n];
            }
            Bit();
        }

   
        public void input_customer_intersection_relationship()
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            strConnection += @"Data Source=|DataDirectory|\Instance_file_Mod.mdb";

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Instance_file_Intersection";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();
            while (myReader001.Read())
            {
                int ddd_true_j = myReader001.GetInt16(1);
                pickupi_match_j_no[ddd_true_j, 1] = myReader001.GetInt16(5);
                dropoffi_match_j_no[ddd_true_j, 1] = myReader001.GetInt16(6);
                pickupi_match_j_no[ddd_true_j, 2] = myReader001.GetInt16(7);
                pickupi_match_j_no[ddd_true_j, 3] = myReader001.GetInt16(8);
                pickupi_match_j_no[ddd_true_j, 4] = myReader001.GetInt16(9);
                pickupi_match_j_no[ddd_true_j, 5] = myReader001.GetInt16(10);
                pickupi_match_j_no[ddd_true_j, 6] = myReader001.GetInt16(11);
                pickupi_match_j_no[ddd_true_j, 7] = myReader001.GetInt16(12);
                pickupi_match_j_no[ddd_true_j, 8] = myReader001.GetInt16(13);
                pickupi_match_j_no[ddd_true_j, 9] = myReader001.GetInt16(14);
                pickupi_match_j_no[ddd_true_j, 10] = myReader001.GetInt16(15);
                dropoffi_match_j_no[ddd_true_j, 2] = myReader001.GetInt16(16);
                dropoffi_match_j_no[ddd_true_j, 3] = myReader001.GetInt16(17);
                dropoffi_match_j_no[ddd_true_j, 4] = myReader001.GetInt16(18);
                dropoffi_match_j_no[ddd_true_j, 5] = myReader001.GetInt16(19);
                dropoffi_match_j_no[ddd_true_j, 6] = myReader001.GetInt16(20);
                dropoffi_match_j_no[ddd_true_j, 7] = myReader001.GetInt16(21);
                dropoffi_match_j_no[ddd_true_j, 8] = myReader001.GetInt16(22);
                dropoffi_match_j_no[ddd_true_j, 9] = myReader001.GetInt16(23);
                dropoffi_match_j_no[ddd_true_j, 10] = myReader001.GetInt16(24);
            }
            myReader001.Close();
            objConnection.Close();
        }

        public void input_data_Manhattan_Intersection_Coordinates_Mod_Time_Car()
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            strConnection += @"Data Source=|DataDirectory|\Manhattan_Intersection_Coordinates_Mod_Time_Car.mdb";


            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Manhattan_Intersection_Coordinates";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();
            while (myReader001.Read())
            {
                int i1 = myReader001.GetInt16(1);
                int i2 = myReader001.GetInt16(2);
                Intersection_car_t_ij[i1, i2] = myReader001.GetDouble(7);
            }
            myReader001.Close();
            objConnection.Close();

        }

        public void input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Car()
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            strConnection += @"Data Source=|DataDirectory|\Manhattan_Intersection_Coordinates_Mod_Dis_Car.mdb";


            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Manhattan_Intersection_Coordinates";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();
            while (myReader001.Read())
            {
                int i1 = myReader001.GetInt16(1);
                int i2 = myReader001.GetInt16(2);
                Intersection_car_d_ij[i1, i2] = myReader001.GetDouble(7);
            }
            myReader001.Close();
            objConnection.Close();

        }

        public void input_data_Manhattan_Intersection_Coordinates_Mod_Time_Pedestrian()
        {
            for (int i1 = 1; i1 < Total_I + 1; i1++)
            {
                for (int i2 = 1; i2 < Total_I + 1; i2++)
                {
                    Intersection_ped_t_ij[i1, i2] = 1000000;
                }
            }

            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";

            strConnection += @"Data Source=|DataDirectory|\Manhattan_Intersection_Coordinates_Mod_Time_Pedestrian.mdb";


            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Manhattan_Intersection_Coordinates";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();
            while (myReader001.Read())
            {
                int i1 = myReader001.GetInt16(1);
                int i2 = myReader001.GetInt16(3);
                Intersection_ped_t_ij[i1, i2] = myReader001.GetDouble(8);


            }
            myReader001.Close();
            objConnection.Close();
        }

        public void input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Pedestrian()
        {
            for (int i1 = 1; i1 < Total_I + 1; i1++)
            {
                for (int i2 = 1; i2 < Total_I + 1; i2++)
                {
                    Intersection_ped_d_ij[i1, i2] = 1000000;
                }
            }

            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";

            strConnection += @"Data Source=|DataDirectory|\Manhattan_Intersection_Coordinates_Mod_Dis_Pedestrian.mdb";


            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Manhattan_Intersection_Coordinates";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();
            while (myReader001.Read())
            {
                int i1 = myReader001.GetInt16(1);
                int i2 = myReader001.GetInt16(3);
                Intersection_ped_d_ij[i1, i2] = myReader001.GetDouble(8);
            }
            myReader001.Close();
            objConnection.Close();
        }

        public void input_five_files_for_intersections_overall_data()
        {

            DateTime begin_total;
            DateTime end_total;
            TimeSpan interval_total;

            begin_total = DateTime.Now;

            input_customer_intersection_relationship();

            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("1) input_customer_intersection_relationship DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));
 
            input_data_Manhattan_Intersection_Coordinates_Mod_Time_Car();

            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("2) input_data_Manhattan_Intersection_Coordinates_Mod_Time_Car DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));
 
            input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Car();

            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("3) input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Car DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));
 
            input_data_Manhattan_Intersection_Coordinates_Mod_Time_Pedestrian();

            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("4) input_data_Manhattan_Intersection_Coordinates_Mod_Time_Pedestrian DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));
 
            input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Pedestrian();

            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("5) input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Pedestrian DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));
 
            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("> > input_five_files_for_intersections_overall_data DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));

        }

     
        public void input_combined_ball_Discre0611_path_s_e_fixdropoff_rolling_0202_yes_flexibledropoff2_matching_nwords_changej(int Start_node, double Start_time, int index_start, int ride_sharing)
        {
            int ll = 0;
            prel_l[ll] = -1;
            d_l[ll] = 1;
            Lall[0] = 1;
            n_l[ll] = Start_node;
            T_l[ll] = 0;
            pipipi_l[ll] = 0;
            S_l[ll] = 0;
            Q_l[ll] = 0;

            for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
            {
                preset1_ll[0, iw] = 0;
                preset2_ll[0, iw] = 0;
                preset3_ll[0, iw] = 0;
            }

            for (int jjj = W + 1; jjj < d1 + 1; jjj++)
            {
                if (status_j[jjj] == 55)
                {
                    if (status_match_j[jjj] == 1)
                    {
                        if (T_l[0] + t_ij[Start_node, jjj + n] > Deadline_j[jjj] - Start_time)
                        {                          
                            preset3_ll[0, lwrd_j[jjj]] = preset3_ll[0, lwrd_j[jjj]] + lbit_j[jjj];                    
                        }
                        else
                        {                        
                            double min_early = int.MaxValue;
                            for (int cc = 1; cc < discre_no + 1; cc++)
                            {
                                double check_hh;
                                if (Real_Coordinate_Time_Pedestrian_data(jjj, cc) > Real_Coordinate_Time_Car_data(Start_node, 1, jjj, cc))
                                {
                                    check_hh = Real_Coordinate_Time_Pedestrian_data(jjj, cc);
                                }
                                else
                                {
                                    check_hh = Real_Coordinate_Time_Car_data(Start_node, 1, jjj, cc);
                                }

                                double HHH = check_hh + (Real_Coordinate_Time_Car_data(jjj, cc, jjj + n, 1));

                                if (HHH < min_early)
                                {
                                    min_early = HHH;
                                }
                            }
                            if (min_early > Deadline_j[jjj] - Start_time)
                            {
                                preset3_ll[0, lwrd_j[jjj]] = preset3_ll[0, lwrd_j[jjj]] + lbit_j[jjj];
                            }                          
                        }
                    }
                }
            }

            if (Start_node != 0 && (preset3_ll[0, lwrd_j[Start_node - n]] & lbit_j[Start_node - n]) != lbit_j[Start_node - n])
            {
                preset3_ll[0, lwrd_j[Start_node - n]] = preset3_ll[0, lwrd_j[Start_node - n]] + lbit_j[Start_node - n];
            }

            int compare_label = 0;

            int oo = 1;
            for (int e = 0; e < Lall[0]; e++)
            {
                if (e == 50000 * oo)
                {
                    Console.WriteLine("Step 1-------  e = " + e + " ll = " + ll + " S_l[ll] = " + S_l[ll]);
                    writer.WriteLine("Step 1-------  e = " + e + " ll = " + ll + " S_l[ll] = " + S_l[ll]);

                    oo = oo + 1;
                }

                if ((d_l[e] == 1 && ride_sharing == 1) || (d_l[e] == 1 && ride_sharing == 0 && S_l[e] <= 2))
                {
                    int i = n_l[e];

                    if (S_l[e] == 0)
                    {
                        for (int j = W + 1; j < d1 + 1; j++) 
                        {
                            if (status_j[j] == 55 || status_j[j] == 88)
                            {
                                if (status_match_j[j] == 1)
                                {
                                    if ((preset3_ll[e, lwrd_j[j]] & lbit_j[j]) != lbit_j[j])
                                    {
                                        ll = ll + 1;
                                        n_l[ll] = j;
                                        T_l[ll] = T_l[e];
                                        prel_l[ll] = e;
                                        Lall[0] = Lall[0] + 1;
                                        d_l[ll] = 1;

                                        for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                        {
                                            preset1_ll[ll, iw] = preset1_ll[e, iw];
                                            preset2_ll[ll, iw] = preset2_ll[e, iw];
                                            preset3_ll[ll, iw] = preset3_ll[e, iw];
                                        }

                                     
                                        preset1_ll[ll, lwrd_j[j]] = preset1_ll[e, lwrd_j[j]] + lbit_j[j];
                                        preset2_ll[ll, lwrd_j[j]] = preset2_ll[e, lwrd_j[j]] + lbit_j[j];
                                        Q_l[ll] = Q_l[e] + q_j[j];
                                        S_l[ll] = S_l[e] + 1;
                                       
                                        for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                        {
                                            preset3_ll[ll, iw] = preset3_ll[e, iw];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else if (S_l[e] != 0)
                    {
                        if ((Q_l[e] == 0 && ride_sharing == 1) || (Q_l[e] == 0 && ride_sharing == 0 && S_l[e] == 2))
                        {
                            d_l[e] = 2;

                            compare_label = compare_label + 1;

                            r = r + 1;
                            for (int jjjj = 1; jjjj < 2 * n + 1; jjjj++)
                            {
                                a_rj[r, jjjj] = 0;
                                seq_rj[r, jjjj] = 0;
                            }
                            start_r[r] = index_start;
                            seq_rj[r, Start_node] = 0 + 999;
                            T_rj[r, Start_node] = 0;
                            end_r[r] = -1;
                            for (int ss = 1; ss < s + 1; ss++)
                            {
                                if (n_location_s[ss] == n_l[e])
                                {
                                    if (time_s[ss] == Start_time / unit + Math.Floor(T_l[e] / unit) + 1)
                                    {
                                        end_r[r] = ss;
                                       
                                        g_r[r] = -pipipi_l[e];
                                        int pre_label = e;
                                        int pre_vertex = n_l[pre_label];
                                        
                                        record_r[r] = S_l[pre_label];            
                                        seq_rj[r, pre_vertex] = (S_l[pre_label] + 999);
                                        T_rj[r, pre_vertex] = (float)T_l[pre_label];
                                        pipipi_rj[r, pre_vertex] = (float)pipipi_l[pre_label];

                                        for (int ii = 0; ii < S_l[e]; ii++)
                                        {
                                            pre_label = prel_l[pre_label];
                                            pre_vertex = n_l[pre_label];
                                            seq_rj[r, pre_vertex] = (S_l[pre_label] + 999);
                                            T_rj[r, pre_vertex] = (float)T_l[pre_label];
                                            pipipi_rj[r, pre_vertex] = (float)pipipi_l[pre_label];

                                            if (pre_vertex <= n && pre_vertex > 0)
                                            {
                                                a_rj[r, pre_vertex] = 1;
                                            }
                                        }
                            
                                     
                                        r_dummy = r_dummy + 1;                           
                                        start_r_dummy[r_dummy] = ss;

                                        for (int ssss = 1; ssss < s + 1; ssss++)
                                        {
                                            if (n_location_s[ssss] == 2 * n + 1)
                                            {
                                                end_r_dummy[r_dummy] = ssss;
                                                break;
                                            }
                                        }
                                        g_r_dummy[r_dummy] = 0;

                                   

                                        break;
                                    }
                                }
                            }
                            if (end_r[r] == -1)
                            {

                                Console.WriteLine("something wrong!!!!!!!!!");
                                Console.WriteLine("Start_node = " + Start_node);
                                Console.WriteLine("start_time = " + Start_time / unit + "end_time = " + T_0 / unit);
                                Console.WriteLine("earliest_time_j[" + (n_l[e] - n) + "] = " + earliest_time_j[n_l[e] - n]);
                                Console.WriteLine("latest_time_j[" + (n_l[e] - n) + "] = " + latest_time_j[n_l[e] - n]);
                                Console.WriteLine("T_l[e] / unit = " + T_l[e] / unit);


                                int pre_label1 = e;
                                int pre_vertex1 = n_l[pre_label1];
                                Console.WriteLine("the last vertex is" + pre_vertex1);
                                for (int ii = 0; ii < S_l[e] - 1; ii++)
                                {
                                    pre_label1 = prel_l[pre_label1];
                                 
                                    pre_vertex1 = n_l[pre_label1];
                                    Console.WriteLine("the last vertex is" + pre_vertex1);
                                
                                }
                            }
                        }
                        

                        if (n_l[e] > n && d_l[e] == 1)
                        {
                            int remark_ll = ll + 1;
                            for (int j = (n + W) + 1; j < (n + d1) + 1; j++)
                            {                               
                                if ((preset1_ll[e, lwrd_j[j - n]] & lbit_j[j - n]) == lbit_j[j - n])
                                {
                                    {
                                        ll = ll + 1;
                                        n_l[ll] = j;
                                        prel_l[ll] = e;
                                        Lall[0] = Lall[0] + 1;
                                        d_l[ll] = 1;
                                        S_l[ll] = S_l[e] + 1;                                      
                                        Q_l[ll] = Q_l[e] + q_j[j];

                                        for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                        {
                                            preset1_ll[ll, iw] = preset1_ll[e, iw];
                                            preset2_ll[ll, iw] = preset2_ll[e, iw];
                                            preset3_ll[ll, iw] = preset3_ll[e, iw];
                                        }
                                        preset1_ll[ll, lwrd_j[j - n]] = preset1_ll[e, lwrd_j[j - n]] - lbit_j[j - n];

                                        One_route_updating_combined_ball_Discre0611_flexibledropoff_changej(ll, Start_time);

                                        if (feasibility == -1 || T_l[ll] > Deadline_j[j - n] - Start_time)
                                        {
                                            d_l[e] = 0; 
                                            ll = remark_ll - 1;
                                            Lall[0] = remark_ll;
                                            break;
                                        }
                                        else
                                        {
                                            if (Q_l[ll] == 0 && d_l[ll] == 1)
                                            {
                                                Dominance_rule_nwords1(ll);
                                            }
                                        }
                                    }
                                }
                            }
                            if (d_l[e] == 1)
                            {
                                for (int j = W + 1; j < d1 + 1; j++)
                                {
                                    if (status_j[j] == 55 || status_j[j] == 88)
                                    {
                                        if (status_match_j[j] == 1)
                                        {
                                            if ((preset3_ll[e, lwrd_j[j]] & lbit_j[j]) != lbit_j[j] && (preset2_ll[e, lwrd_j[j]] & lbit_j[j]) != lbit_j[j] && Q_l[e] + q_j[j] <= Q)
                                            {
                                                ll = ll + 1;
                                                n_l[ll] = j;                                     
                                                T_l[ll] = T_l[e];
                                                prel_l[ll] = e;
                                                Lall[0] = Lall[0] + 1;
                                                d_l[ll] = 1;
                                                S_l[ll] = S_l[e] + 1;
                                                Q_l[ll] = Q_l[e] + q_j[j];

                                                for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                                {
                                                    preset1_ll[ll, iw] = preset1_ll[e, iw];
                                                    preset2_ll[ll, iw] = preset2_ll[e, iw];
                                                    preset3_ll[ll, iw] = preset3_ll[e, iw];
                                                }
                                                preset1_ll[ll, lwrd_j[j]] = preset1_ll[e, lwrd_j[j]] + lbit_j[j];
                                                preset2_ll[ll, lwrd_j[j]] = preset2_ll[e, lwrd_j[j]] + lbit_j[j];
                                            }
                                        }
                                    }

                                }

                            }
                        }
                        else if (n_l[e] <= n && d_l[e] == 1)
                        {
                           
                            int remark_ll = ll + 1;
                            for (int j = (n + W) + 1; j < (n + d1) + 1; j++)
                            {
                                if ((preset1_ll[e, lwrd_j[j - n]] & lbit_j[j - n]) == lbit_j[j - n])
                                {
                                    ll = ll + 1;
                                    n_l[ll] = j;
                                    prel_l[ll] = e;
                                    Lall[0] = Lall[0] + 1;
                                    d_l[ll] = 1;
                                    S_l[ll] = S_l[e] + 1;
                                    Q_l[ll] = Q_l[e] + q_j[j];

                                    for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                    {
                                        preset1_ll[ll, iw] = preset1_ll[e, iw];
                                        preset2_ll[ll, iw] = preset2_ll[e, iw];
                                        preset3_ll[ll, iw] = preset3_ll[e, iw];
                                    }
                                    preset1_ll[ll, lwrd_j[j - n]] = preset1_ll[ll, lwrd_j[j - n]] - lbit_j[j - n];

                                    One_route_updating_combined_ball_Discre0611_flexibledropoff_changej(ll, Start_time);

                                    if (feasibility == -1 || T_l[ll] > Deadline_j[j - n] - Start_time)
                                    {
                                        d_l[ll] = 0;
                                        ll = remark_ll - 1;
                                        Lall[0] = remark_ll;
                                        d_l[e] = 0; 
                                        for (int llll = e + 1; llll < Lall[0] + 1; llll++)
                                        {
                                            if (d_l[llll] == 1)
                                            {
                                                if (prel_l[llll] == prel_l[e])
                                                {
                                                    if ((preset3_ll[llll, lwrd_j[n_l[e]]] & lbit_j[n_l[e]]) != lbit_j[n_l[e]])
                                                    {
                                                        preset3_ll[llll, lwrd_j[n_l[e]]] = preset3_ll[llll, lwrd_j[n_l[e]]] + lbit_j[n_l[e]];
                                                    }
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    else
                                    {
                                        if (Q_l[ll] == 0 && d_l[ll] == 1)
                                        {
                                            Dominance_rule_nwords1(ll);
                                        }
                                    }
                                }
                            }

                            if (d_l[e] == 1)
                            {
                                for (int j = W + 1; j < d1 + 1; j++)
                                {
                                    if (status_j[j] == 55 || status_j[j] == 88)
                                    {
                                        if (status_match_j[j] == 1)
                                        {
                                            if ((preset3_ll[e, lwrd_j[j]] & lbit_j[j]) != lbit_j[j] && (preset2_ll[e, lwrd_j[j]] & lbit_j[j]) != lbit_j[j] && Q_l[e] + q_j[j] <= Q)
                                            {
                                                ll = ll + 1;
                                                n_l[ll] = j;                                              
                                                T_l[ll] = T_l[e];
                                                prel_l[ll] = e;
                                                Lall[0] = Lall[0] + 1;
                                                d_l[ll] = 1;
                                                S_l[ll] = S_l[e] + 1;
                                                Q_l[ll] = Q_l[e] + q_j[j];

                                                for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                                {
                                                    preset1_ll[ll, iw] = preset1_ll[e, iw];
                                                    preset2_ll[ll, iw] = preset2_ll[e, iw];
                                                    preset3_ll[ll, iw] = preset3_ll[e, iw];
                                                }
                                                preset1_ll[ll, lwrd_j[j]] = preset1_ll[ll, lwrd_j[j]] + lbit_j[j];
                                                preset2_ll[ll, lwrd_j[j]] = preset2_ll[ll, lwrd_j[j]] + lbit_j[j];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        public double eps = 1e-1;

        public void output_information_route_Discre0611_20210104_rr_rolling0202_ball_flexibledropoff_clearversion(int rr)
        {

          
            eps = 0.001;
            mm_eps = 0.001;
            check_opt = 1;

        
            int last_pre_vertex = -1000;
            int last_last_pre_vertex = -1000;

            for (int ii = record_r[rr]; ii > 0; ii--)
            {
                 last_pre_vertex = n_location_s[end_r[rr]];
                for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++) 
                {
                    if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - 1)
                    {
                        last_last_pre_vertex = jjjj;
                        break;
                    }
                }
                 
                int bb = 1;
                while (bb < ii)
                {
                    for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                    {
                        if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - bb)
                        {
                            last_pre_vertex = jjjj;
                            break;
                        }
                    }
                    bb = bb + 1;

                    for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                    {
                        if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - bb)
                        {
                            last_last_pre_vertex = jjjj;
                            break;
                        }
                    }
                }

            }

            One_route_updating_combined_ball_Discre0611_20210104_rr_flexibledropoff(rr, last_pre_vertex);

         

            double[] determined_v_i = new double[2 * nprime + 2];

            last_pre_vertex = n_location_s[end_r[rr]];

            for (int ii = 0; ii < record_r[rr] - 1; ii++)
            {
               

                for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                {
                    if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - ii - 1)
                    {
                        last_pre_vertex = jjjj;                      
                        break;
                    }
                }           
            }


            determined_v_i[0] = time_s[start_r[rr]] * unit;

            if (n_location_s[start_r[rr]] == 0)
            {
                aa_ln[0] = aao_j[0];
                bb_ln[0] = bbo_j[0];

                aabb_ln_data[0] = 0;
                aabb_ln_data_candidate_no[0] = 1;
              
            }
            else
            {
                aa_ln[0] = aad_j[n_location_s[start_r[rr]] - n];
                bb_ln[0] = bbd_j[n_location_s[start_r[rr]] - n];

                aabb_ln_data[0] = n_location_s[start_r[rr]];
                aabb_ln_data_candidate_no[0] = 1;
               
            }


           


            for (int ii = 1; ii < record_r[rr] + 1; ii++)
            {
                int bb = 1;
                last_pre_vertex = n_location_s[end_r[rr]];
                while (bb <= record_r[rr] - ii)
                {

                    for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                    {
                        if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - bb)
                        {
                            last_pre_vertex = jjjj;
                            break;
                        }
                    }
                    bb = bb + 1;
                }

               
                if (last_pre_vertex <= n && last_pre_vertex > 0)
                {
                    if (determined_v_i[0] + Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]) > determined_v_i[ii - 1] + Real_Coordinate_Time_Car_data(aabb_ln_data[ii - 1], aabb_ln_data_candidate_no[ii - 1], aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]))
                    {
                        determined_v_i[ii] = determined_v_i[0] + Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);
                    }
                    else
                    {
                        determined_v_i[ii] = determined_v_i[ii - 1] + Real_Coordinate_Time_Car_data(aabb_ln_data[ii - 1], aabb_ln_data_candidate_no[ii - 1], aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);
                    }
                }
                else
                {
                    determined_v_i[ii] = determined_v_i[ii - 1] + Real_Coordinate_Time_Car_data(aabb_ln_data[ii - 1], aabb_ln_data_candidate_no[ii - 1], aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);
                }

            }

            eps = 1e-1;
            mm_eps = 0.1;
        }

     
      
        public int vehicle_stop = 0;
        public double vehicle_profit = 0;
    


     
        public void output_information_route_Discre0611_20210104_rr_to_excel_flexibledropoff_rolling0202_0212_ball_matching_clearversion(int vehicle, int sub_route, int rr, double next_resolve_end)
        {

            if (time_s[start_r[rr]] * unit < next_resolve_end)
            {
                sub_path = sub_path + 1;
                vehicle_subpath[sub_path] = vehicle;

                time_start_subpath[sub_path] = time_s[start_r[rr]] * unit;
                n_location_start_subpath[sub_path] = n_location_s[start_r[rr]];

                time_end_subpath[sub_path] = time_s[end_r[rr]] * unit; 
                n_location_end_subpath[sub_path] = n_location_s[end_r[rr]];

                profit_subpath[sub_path] = g_r[rr];

             
                if (set_s[start_r[rr]] == -1)
                {
                    vehicle_profit = 0;
                    vehicle_stop = 0;
                }

                double[] determined_v_i = new double[2 * nprime + 2];

                determined_v_i[0] = time_s[start_r[rr]] * unit;

                int last_pre_vertex;
                for (int ii = 1; ii < record_r[rr] + 1; ii++)
                {
                    int bb = 1;
                    last_pre_vertex = n_location_s[end_r[rr]];
                    while (bb <= record_r[rr] - ii)
                    {
                        for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                        {
                            if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - bb) 
                            {
                                last_pre_vertex = jjjj;
                                break;
                            }
                        }
                        bb = bb + 1;
                    }

                    if (last_pre_vertex <= n && last_pre_vertex > 0)
                    {
                        if (determined_v_i[0] + Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]) > determined_v_i[ii - 1] + Real_Coordinate_Time_Car_data(aabb_ln_data[ii - 1], aabb_ln_data_candidate_no[ii - 1], aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]))
                        {
                            determined_v_i[ii] = determined_v_i[0] + Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);
                        }
                        else
                        {
                            determined_v_i[ii] = determined_v_i[ii - 1] + Real_Coordinate_Time_Car_data(aabb_ln_data[ii - 1], aabb_ln_data_candidate_no[ii - 1], aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);
                        }

                        vehicle_stop = vehicle_stop + 1;

                        vehicle_profit = vehicle_profit + g_j[last_pre_vertex];

                        a_subpathj[sub_path, last_pre_vertex] = 1;

                  
                        stop_vehicle[vehicle] = stop_vehicle[vehicle] + 1;

                        truej_vehicle_stop[vehicle, stop_vehicle[vehicle]] = true_j[last_pre_vertex];
                        subpath_vehicle_stop[vehicle, stop_vehicle[vehicle]] = sub_path;

                        release_vehicle_truej[true_j[last_pre_vertex]] = vehicle;
                        release_stop_seq_truej[true_j[last_pre_vertex]] = stop_vehicle[vehicle];

                        release_subpath_truej[true_j[last_pre_vertex]] = sub_path;

                        q_truej[true_j[last_pre_vertex]] = q_j[last_pre_vertex];
                        g_truej[true_j[last_pre_vertex]] = g_j[last_pre_vertex];

                        aao_truej[true_j[last_pre_vertex]] = aao_j[last_pre_vertex];
                        bbo_truej[true_j[last_pre_vertex]] = bbo_j[last_pre_vertex];

               

                        release_aao_truej[true_j[last_pre_vertex]] = aac_i[pickupi_match_j_no[true_j[aabb_ln_data[ii]], aabb_ln_data_candidate_no[ii]]];
                        release_bbo_truej[true_j[last_pre_vertex]] = bbc_i[pickupi_match_j_no[true_j[aabb_ln_data[ii]], aabb_ln_data_candidate_no[ii]]];

                        release_aabbo_truej_data_candidate_no[true_j[last_pre_vertex]] = aabb_ln_data_candidate_no[ii];

                        release_start_time_truej[true_j[last_pre_vertex]] = time_s[start_r[rr]] * unit;
                        release_walk_time_truej[true_j[last_pre_vertex]] = (determined_v_i[0] + Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]));

                        release_arrivestop_time_truej[true_j[last_pre_vertex]] = determined_v_i[ii];
                        release_arrivehome_time_truej[true_j[last_pre_vertex]] = -1;

                    }
                    else
                    {
                        determined_v_i[ii] = determined_v_i[ii - 1] + Real_Coordinate_Time_Car_data(aabb_ln_data[ii - 1], aabb_ln_data_candidate_no[ii - 1], aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);


                        vehicle_stop = vehicle_stop + 1;

           
                        stop_vehicle[vehicle] = stop_vehicle[vehicle] + 1;

                        truej_vehicle_stop[vehicle, stop_vehicle[vehicle]] = true_j[last_pre_vertex - n] + truen;
                        subpath_vehicle_stop[vehicle, stop_vehicle[vehicle]] = sub_path;

                        release_vehicle_truej[true_j[last_pre_vertex - n] + truen] = vehicle;
                        release_stop_seq_truej[true_j[last_pre_vertex - n] + truen] = stop_vehicle[vehicle];

                        release_subpath_truej[true_j[last_pre_vertex - n] + truen] = sub_path;

                        q_truej[true_j[last_pre_vertex - n] + truen] = q_j[last_pre_vertex];

                        aad_truej[true_j[last_pre_vertex - n]] = aad_j[last_pre_vertex - n];
                        bbd_truej[true_j[last_pre_vertex - n]] = bbd_j[last_pre_vertex - n];

                     
                        release_aad_truej[true_j[last_pre_vertex - n]] = aac_i[dropoffi_match_j_no[true_j[aabb_ln_data[ii] - n], aabb_ln_data_candidate_no[ii]]];
                        release_bbd_truej[true_j[last_pre_vertex - n]] = bbc_i[dropoffi_match_j_no[true_j[aabb_ln_data[ii] - n], aabb_ln_data_candidate_no[ii]]];

                        release_aabbd_truej_data_candidate_no[true_j[last_pre_vertex - n]] = aabb_ln_data_candidate_no[ii];


                        release_start_time_truej[true_j[last_pre_vertex - n] + truen] = time_s[start_r[rr]] * unit;
                        release_walk_time_truej[true_j[last_pre_vertex - n] + truen] = Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]);

                        release_arrivestop_time_truej[true_j[last_pre_vertex - n] + truen] = determined_v_i[ii];
                        release_arrivehome_time_truej[true_j[last_pre_vertex - n] + truen] = (determined_v_i[ii] + Real_Coordinate_Time_Pedestrian_data(aabb_ln_data[ii], aabb_ln_data_candidate_no[ii]));


                        deadline_truej[true_j[last_pre_vertex - n]] = Deadline_j[last_pre_vertex - n];

                    }
                }

                 fare_subpath[sub_path] = vehicle_profit;
 
                vehicle_profit = vehicle_profit - (determined_v_i[record_r[rr]] - determined_v_i[0]) * cost;
 
                detail_time_end_subpath[sub_path] = determined_v_i[record_r[rr]];
                detail_duration_subpath[sub_path] = (determined_v_i[record_r[rr]] - determined_v_i[0]);
                cost_subpath[sub_path] = (determined_v_i[record_r[rr]] - determined_v_i[0]) * cost;
 
                for (int rrr = 1; rrr < r + 1; rrr++)
                {
                    if (par_z_r[rrr] == 1 && n_location_s[start_r[rrr]] == n_location_s[end_r[rr]] && n_location_s[end_r[rrr]] == (2 * n + 1))
                    {
                        vehicle_stop = vehicle_stop + 1;

                        sub_route = sub_route + 1;

                        break;
                    }
                }
            }
            else
            {
                int last_pre_vertex;
                for (int ii = 1; ii < record_r[rr] + 1; ii++)
                {
                    int bb = 1;
                    last_pre_vertex = n_location_s[end_r[rr]];
                    while (bb <= record_r[rr] - ii)
                    {
                        for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                        {
                            if (seq_rj[rr, jjjj] == seq_rj[rr, n_location_s[end_r[rr]]] - bb)
                            {
                                last_pre_vertex = jjjj;
                                break;
                            }
                        }
                        bb = bb + 1;
                    }

                    if (last_pre_vertex <= n && last_pre_vertex > 0)
                    {
                        cust_set_j[last_pre_vertex] = vehicle;
                        if (first_cust_v[vehicle] == 0)
                        {
                            first_cust_v[vehicle] = last_pre_vertex;
                        }
                    }
                }
            }
        }

      
   
        public const int Com_cust_set_number = 2;
        public int[,] Cust_set_j1_number = new int[n + 1, Com_cust_set_number + 1];
        public int[,] Rank_vvv1_number = new int[W + 1, W + 1];

        public int[] status_match_j = new int[n + 1];

     
        public int KKK = 50;

        public void create_time_space_within_timewindow_for_rolling_0202_ball_flexibledropoff_matching_changej_0405(int KKK, double resolve_end, int no_rolling)
        {

            int NOO = 500;

            double[,] travel_vvv1vvv2 = new double[W + 1, W + 1];
            double[,] travel_vvvjjj = new double[W + 1, n + 1];

            for (int vvv1 = 1; vvv1 < W + 1; vvv1++)
            {
                travel_vvv1vvv2[vvv1, vvv1] = 0;
                for (int vvv2 = vvv1 + 1; vvv2 < W + 1; vvv2++)
                {
                    travel_vvv1vvv2[vvv1, vvv2] = Math.Abs(EarliestTime_v[vvv1] - EarliestTime_v[vvv2]) + (Real_Coordinate_Time_Car_data(n_location_v[vvv1], 1, n_location_v[vvv2], 1) + Real_Coordinate_Time_Car_data(n_location_v[vvv2], 1, n_location_v[vvv1], 1)) / 2.0;
                    travel_vvv1vvv2[vvv2, vvv1] = travel_vvv1vvv2[vvv1, vvv2];
                }
            }

            for (int vvv1 = 1; vvv1 < W + 1; vvv1++)
            {
                double previous_min = int.MinValue;

                int[] occupy = new int[W + 1];

                Rank_vvv1_number[vvv1, 1] = vvv1;
                occupy[vvv1] = 1;

                for (int no = 2; no < NOO + 1; no++)
                {
                    double min = int.MaxValue;
                    int min_vvv2 = -111;

                    for (int vvv2 = 1; vvv2 < W + 1; vvv2++)
                    {
                        if (min > travel_vvv1vvv2[vvv1, vvv2] && travel_vvv1vvv2[vvv1, vvv2] >= previous_min && occupy[vvv2] == 0)
                        {
                            min = travel_vvv1vvv2[vvv1, vvv2];
                            min_vvv2 = vvv2;
                        }
                    }
                     previous_min = min;
                    Rank_vvv1_number[vvv1, no] = min_vvv2;
                    occupy[min_vvv2] = 1;
                 }
            }


 
            int[] minvvv_jjj = new int[n + 1];
            for (int jjj = W + 1; jjj < d1 + 1; jjj++) 
            {
                minvvv_jjj[jjj] = 0; 
                if (status_j[jjj] == 55)
                {
                    double min = int.MaxValue;
                    int min_vvv = -111;

                    for (int vvv = 1; vvv < W + 1; vvv++)
                    {
                        travel_vvvjjj[vvv, jjj] = EarliestTime_v[vvv] + (Real_Coordinate_Time_Car_data(n_location_v[vvv], 1, jjj, 1) + Real_Coordinate_Time_Car_data(jjj, 1, n_location_v[vvv], 1)) / 2.0;
                        if (min > travel_vvvjjj[vvv, jjj])
                        {
                            min = travel_vvvjjj[vvv, jjj];
                            min_vvv = vvv;
                        
                        }
                    }
                    minvvv_jjj[jjj] = min_vvv;
                 }
            }
        
            for (int w = 1; w < W + 1; w++) 
            {
                s = s + 1;
                check_s[s] = 1;
                n_location_s[s] = n_location_v[w];
                time_s[s] = (int)(EarliestTime_v[w] / unit);
               
                set_s[s] = -1;
                vehicle_s[s] = w;
 
            }
 
            {
                s = s + 1;
                check_s[s] = 1;
                n_location_s[s] = 2 * n + 1;
                time_s[s] = (int)Math.Floor(T_0 / unit);
                set_s[s] = 1;
            }


            for (int j = W + 1; j < d1 + 1; j++)
            {
                if (status_j[j] == 55 || status_j[j] == 88)
                {
                   
                    double min_early = int.MaxValue;
                    for (int w = 1; w < W + 1; w++)
                    {
                        for (int cc = 1; cc < discre_no + 1; cc++)
                        {
                            double check_hh;
                            if (Real_Coordinate_Time_Pedestrian_data(j, cc) > Real_Coordinate_Time_Car_data(n_location_v[w], 1, j, cc))
                            {
                                check_hh = Real_Coordinate_Time_Pedestrian_data(j, cc);
                            }
                            else
                            {
                                check_hh = Real_Coordinate_Time_Car_data(n_location_v[w], 1, j, cc);
                            }

                            double HHH = (int)(EarliestTime_v[w] / unit) * unit + check_hh + Real_Coordinate_Time_Car_data(j, cc, j + n, 1);

                            if (HHH < min_early)
                            {
                                min_early = HHH;
                            }
                        }
                    }

                    earliest_time_j[j] = (int)Math.Floor(min_early / unit) - 10; 

                    latest_time_j[j] = (int)Math.Ceiling(Deadline_j[j] / unit);

                    if (latest_time_j[j] >= earliest_time_j[j])
                    {
                        for (int t = earliest_time_j[j]; t < latest_time_j[j] + 1; t++)
                        {
                            s = s + 1;

                            if (t >= earliest_time_j[j] + 10)
                            {
                                check_s[s] = 1;
                            }
                            else
                            {
                                check_s[s] = -1;
                            }

                            n_location_s[s] = j + n;
                            time_s[s] = t;

                         
                        }
                    }
                    else
                    {
                        status_j[j] = -88;
                        Console.WriteLine("true_j = " + true_j[j] + " >> customer " + j + " is rejected, because latest_time_j (" + latest_time_j[j] + ") < earliest_time_j (" + earliest_time_j[j] + ")");
                        direct_rejected_new_customer_roll[no_rolling] = direct_rejected_new_customer_roll[no_rolling] + 1;
                    }
                }

        
            }

            Console.WriteLine(" direct_rejected_new_customer_roll[" + no_rolling + " ] = " + direct_rejected_new_customer_roll[no_rolling]);

            Console.WriteLine("In total:: :: ::: ::: s = " + s);
            writer.WriteLine("In total:: :: ::: ::: s = " + s);

            timespace_node_roll[no_rolling] = s;
            Console.WriteLine(" timespace_node_roll[" + no_rolling + "] = " + timespace_node_roll[no_rolling]);

            t_ij[0, 2 * n + 1] = 0;
 
            int Time_longest_dp = -1;
            int Num_dp_longest = -1;

            DateTime begin_dp;
            DateTime end_dp;
            TimeSpan interval_dp;

            DateTime begin_dp_total;
            DateTime end_dp_total;
            TimeSpan interval_dp_total;

            begin_dp_total = DateTime.Now;

            int num_dp = 0;
            for (int s_1 = 1; s_1 < s + 1; s_1++)
            {
                if (n_location_s[s_1] != 2 * n + 1 && check_s[s_1] == 1) 
                {
                    num_dp = num_dp + 1;
  

                    begin_dp = DateTime.Now;

                    int accum_kkk = 0;
 

                    for (int j = W + 1; j < d1 + 1; j++)
                    {
                        status_match_j[j] = 0;
                        if (status_j[j] == 55)
                        {
                            status_match_j[j] = 1;
                            accum_kkk = accum_kkk + 1;
                        }
                    }

                    if (set_s[s_1] == -1)
                    {
                        int vvvv1 = vehicle_s[s_1];


                        for (int xx = 1; xx < NOO + 1; xx++)
                        {
                            for (int jjj = W + 1; jjj < d1 + 1; jjj++)
                            {
                                if (status_j[jjj] == 88 && cust_set_j[jjj] == Rank_vvv1_number[vvvv1, xx])
                                {
                                    status_match_j[jjj] = 2;
                                    accum_kkk = accum_kkk + 1;
                                }
                            }

                            if (xx == 1 || accum_kkk < KKK)
                            {
                                for (int jjj = W + 1; jjj < d1 + 1; jjj++)
                                {
                                    if (status_match_j[jjj] == 2)
                                    {
                                        status_match_j[jjj] = 1;
                                     }
                                }
                            }
                            else
                            {
                                break;
                            }
                        }

                    }
                    else 
                    {
                        int j1 = n_location_s[s_1] - n;

                        int vvvv1 = 0;

                        if (status_j[j1] == 88)
                        {
                            vvvv1 = cust_set_j[j1];
                        }
                        if (status_j[j1] == 55)
                        {
                            vvvv1 = minvvv_jjj[j1];
                        }


                        for (int xx = 1; xx < NOO + 1; xx++)
                        {
                            for (int jjj = W + 1; jjj < d1 + 1; jjj++)
                            {
                                if (status_j[jjj] == 88 && cust_set_j[jjj] == Rank_vvv1_number[vvvv1, xx])
                                {
                                    status_match_j[jjj] = 2;
                                    accum_kkk = accum_kkk + 1;
                                }
                            }

                            if (xx == 1 || accum_kkk < KKK)
                            {
                                for (int jjj = W + 1; jjj < d1 + 1; jjj++)
                                {
                                    if (status_match_j[jjj] == 2)
                                    {
                                        status_match_j[jjj] = 1;
                                    }
                                }
                            }
                            else
                            {
                                break;
                            }
                        }

                    }


                  
                    input_combined_ball_Discre0611_path_s_e_fixdropoff_rolling_0202_yes_flexibledropoff2_matching_nwords_changej(n_location_s[s_1], time_s[s_1] * unit, s_1, ride_sharing);


                    end_dp = DateTime.Now;
                    interval_dp = end_dp.Subtract(begin_dp);
                   

                    if (Convert.ToInt32(interval_dp.TotalMilliseconds) > Time_longest_dp)
                    {
                        Time_longest_dp = Convert.ToInt32(interval_dp.TotalMilliseconds);
                        Num_dp_longest = num_dp;
                    }
                }
            }
            Console.WriteLine("num_dp = " + num_dp);

            timespace_node_roll[no_rolling] = num_dp;
            Console.WriteLine("Update : timespace_node_roll[" + no_rolling + "] = " + timespace_node_roll[no_rolling]);

            end_dp_total = DateTime.Now;
            interval_dp_total = end_dp_total.Subtract(begin_dp_total);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("TotalTotalTotalTotal >> DP_total DONE! The Totmin is " + Convert.ToInt32(interval_dp_total.TotalMinutes) + " or Totsec: " + Convert.ToInt32(interval_dp_total.TotalSeconds));
            writer.WriteLine("*******************************************************************");
            writer.WriteLine("TotalTotalTotalTotal >> DP_total DONE! The Totmin is " + Convert.ToInt32(interval_dp_total.TotalMinutes) + " or Totsec: " + Convert.ToInt32(interval_dp_total.TotalSeconds));

            Console.WriteLine("*******************************************************************");
            Console.WriteLine("Parallelize >> DP_total DONE! The Totsec is " + (int)(Time_longest_dp / 1000.0) + " or TotalMilliseconds: " + Time_longest_dp);
            writer.WriteLine("*******************************************************************");
            writer.WriteLine("Parallelize >> DP_total DONE! The Totsec is " + (int)(Time_longest_dp / 1000.0) + " or TotalMilliseconds: " + Time_longest_dp);

            DP_CPU_time_roll[no_rolling] = Time_longest_dp;
            Console.WriteLine("DP_CPU_time_roll[" + no_rolling + "] = " + DP_CPU_time_roll[no_rolling]);
           
            
        }

       
   

        public int[] par_z_r = new int[R + 1];
        
   

        public int[] cust_set_j = new int[n + 1];
        public int[] first_cust_v = new int[W + 1];
     
   
    
        public void deterministic_fomulation_Discre_rolling_0202_ball_flexibledropoff_matching_clearversion_changej_improve(int W, double next_resolve_end, int no_rolling)
        {
            DateTime begin_cplex0;
            DateTime end_cplex0;
            TimeSpan interval_cplex0;

            begin_cplex0 = DateTime.Now;

            Console.WriteLine("vehicle number W = " + W);
            writer.WriteLine("vehicle number W = " + W);

            Cplex Multiple_DAR = new Cplex();
 
            Multiple_DAR.SetOut(null);
 
            int max_times_seconds_tununing = 3600;
            Multiple_DAR.SetParam(ILOG.CPLEX.Cplex.DoubleParam.TiLim, max_times_seconds_tununing);
            Multiple_DAR.SetParam(ILOG.CPLEX.Cplex.DoubleParam.EpGap, 0.0001);

            INumVar[] empty_0 = new INumVar[1];
            empty_0 = Multiple_DAR.NumVarArray(1, 0, 0, NumVarType.Float);

            string[] var_z_r = new string[r + 1];
            INumVar[] z_r = new INumVar[r + 1];
            for (int xx = 0; xx < r + 1; xx++)
            {
                var_z_r[xx] = string.Format("z_r({0})", xx);
            }
            z_r = Multiple_DAR.NumVarArray(r + 1, 0, 1, NumVarType.Bool, var_z_r);

            string[] var_z_r_dummy = new string[r_dummy + 1];
            INumVar[] z_r_dummy = new INumVar[r_dummy + 1];
            for (int xx = 0; xx < r_dummy + 1; xx++)
            {
                var_z_r_dummy[xx] = string.Format("z_r_dummy({0})", xx);
            }
            z_r_dummy = Multiple_DAR.NumVarArray(r_dummy + 1, 0, 1, NumVarType.Bool, var_z_r_dummy);

            INumVar[] vehicle_no = new INumVar[1];
            vehicle_no = Multiple_DAR.NumVarArray(1, 0, W, NumVarType.Int);

            INumExpr[] expr2 = new INumExpr[1];
            INumExpr[] expr3 = new INumExpr[1];
            INumExpr[] expr4 = new INumExpr[1];
            INumExpr[] expr5 = new INumExpr[1];
            INumExpr[] expr45 = new INumExpr[1];
            expr45[0] = Multiple_DAR.Prod(0, empty_0[0]);
            INumExpr[] expr6 = new INumExpr[1];
            INumExpr[] expr7 = new INumExpr[1];
            for (int ss = 1; ss < s + 1; ss++)
            {
                 if (set_s[ss] != 1 && set_s[ss] != -1)
                {
                    expr2[0] = Multiple_DAR.Prod(0, empty_0[0]);
                    expr3[0] = Multiple_DAR.Prod(0, empty_0[0]);

                    for (int rr = 1; rr < r + 1; rr++)
                    {
                        if (start_r[rr] == ss && set_s[end_r[rr]] != -1)  
                        {
                            expr2[0] = Multiple_DAR.Sum(expr2[0], z_r[rr]);
                        }
                        if (end_r[rr] == ss && set_s[start_r[rr]] != 1)   
                        {
                            expr3[0] = Multiple_DAR.Sum(expr3[0], z_r[rr]);
                        }
                    }

                    for (int rr = 1; rr < r_dummy + 1; rr++) 
                    {
                        if (start_r_dummy[rr] == ss && set_s[end_r_dummy[rr]] != -1) 
                        {
                            expr2[0] = Multiple_DAR.Sum(expr2[0], z_r_dummy[rr]);
                        }
                        if (end_r_dummy[rr] == ss && set_s[start_r_dummy[rr]] != 1)  
                        {
                            expr3[0] = Multiple_DAR.Sum(expr3[0], z_r_dummy[rr]);
                        }
                    }

                    Multiple_DAR.AddEq(expr2[0], expr3[0]);
                }
                
                if (set_s[ss] == -1)
                {
                     expr4[0] = Multiple_DAR.Prod(0, empty_0[0]);
                    expr5[0] = Multiple_DAR.Prod(0, empty_0[0]);
                    for (int rr = 1; rr < r + 1; rr++)
                    {
                        if (start_r[rr] == ss && set_s[end_r[rr]] != -1) 
                        {
                            expr4[0] = Multiple_DAR.Sum(expr4[0], z_r[rr]);
                        }
                        if (end_r[rr] == ss && set_s[start_r[rr]] != 1) 
                        {
                            expr5[0] = Multiple_DAR.Sum(expr5[0], z_r[rr]);
                        }
                    }
                     Multiple_DAR.AddLe(Multiple_DAR.Sum(expr4[0], Multiple_DAR.Prod(-1, expr5[0])), 1); 
                    expr45[0] = Multiple_DAR.Sum(expr45[0], Multiple_DAR.Sum(expr4[0], Multiple_DAR.Prod(-1, expr5[0]))); 
                }
                 if (set_s[ss] == 1)
                {
 
                    expr6[0] = Multiple_DAR.Prod(0, empty_0[0]);
                    expr7[0] = Multiple_DAR.Prod(0, empty_0[0]);
                     for (int rr = 1; rr < r_dummy + 1; rr++)  
                    {
                        if (start_r_dummy[rr] == ss && set_s[end_r_dummy[rr]] != -1) 
                        {
                            expr6[0] = Multiple_DAR.Sum(expr6[0], z_r_dummy[rr]);
                        }
                        if (end_r_dummy[rr] == ss && set_s[start_r_dummy[rr]] != 1) 
                        {
                            expr7[0] = Multiple_DAR.Sum(expr7[0], z_r_dummy[rr]);
                        }
                    }
                    Multiple_DAR.AddEq(Multiple_DAR.Sum(expr6[0], Multiple_DAR.Prod(-1, expr7[0])), Multiple_DAR.Prod(-1, vehicle_no[0]));
                }
            }

            Multiple_DAR.AddEq(expr45[0], vehicle_no[0]);
            Multiple_DAR.AddLe(vehicle_no[0], W);


            INumExpr[] expr8 = new INumExpr[1];
            INumExpr[] expr8_2 = new INumExpr[1];
            for (int j = W + 1; j < d1 + 1; j++)
            {
                if (status_j[j] == 55)
                {
                    expr8[0] = Multiple_DAR.Prod(0, empty_0[0]);
                    for (int rr = 1; rr < r + 1; rr++)
                    {
                        expr8[0] = Multiple_DAR.Sum(expr8[0], Multiple_DAR.Prod(z_r[rr], a_rj[rr, j]));
                    }
                    Multiple_DAR.AddLe(expr8[0], 1);
                }

                if (status_j[j] == 88)
                {
                     expr8_2[0] = Multiple_DAR.Prod(0, empty_0[0]);
                    for (int rr = 1; rr < r + 1; rr++)
                    {
                        expr8_2[0] = Multiple_DAR.Sum(expr8_2[0], Multiple_DAR.Prod(z_r[rr], a_rj[rr, j]));
                    }
                    Multiple_DAR.AddEq(expr8_2[0], 1);
                }
            }

            string[] var_obj = new string[1];
            INumVar[] obj = new INumVar[1];
            for (int xx = 0; xx < 1; xx++)
            {
                var_obj[xx] = string.Format("obj({0})", xx);
            }
            obj = Multiple_DAR.NumVarArray(1, -int.MaxValue, int.MaxValue, NumVarType.Float, var_obj);


            INumExpr[] expr9 = new INumExpr[1];
            expr9[0] = Multiple_DAR.Prod(0, empty_0[0]);
            for (int rr = 1; rr < r + 1; rr++)
            {
                expr9[0] = Multiple_DAR.Sum(expr9[0], Multiple_DAR.Prod(z_r[rr], g_r[rr]));
            }

            Multiple_DAR.AddEq(expr9[0], obj[0]);

            Multiple_DAR.AddMaximize(obj[0]);


            end_cplex0 = DateTime.Now;
            interval_cplex0 = end_cplex0.Subtract(begin_cplex0);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("CPLEX Read constraint! The Total Time is " + Convert.ToString(interval_cplex0.Hours) + ":" + Convert.ToString(interval_cplex0.Minutes) + ":" + Convert.ToString(interval_cplex0.Seconds));
            writer.WriteLine("*******************************************************************");
            writer.WriteLine("CPLEX Read constraint! The Total Time is " + Convert.ToString(interval_cplex0.Hours) + ":" + Convert.ToString(interval_cplex0.Minutes) + ":" + Convert.ToString(interval_cplex0.Seconds));

            if (Multiple_DAR.Solve())
            {
                end_cplex0 = DateTime.Now;
                interval_cplex0 = end_cplex0.Subtract(begin_cplex0);
                Console.WriteLine("*******************************************************************");
                Console.WriteLine("CPLEX Read + Solve! The Total Time is " + Convert.ToString(interval_cplex0.Hours) + ":" + Convert.ToString(interval_cplex0.Minutes) + ":" + Convert.ToString(interval_cplex0.Seconds));
                writer.WriteLine("*******************************************************************");
                writer.WriteLine("CPLEX Read + Solve! The Total Time is " + Convert.ToString(interval_cplex0.Hours) + ":" + Convert.ToString(interval_cplex0.Minutes) + ":" + Convert.ToString(interval_cplex0.Seconds));

                double Multiple_DAR_benefit = Multiple_DAR.ObjValue;
                Console.WriteLine("$$$The Multiple_DAR_benefit is " + Multiple_DAR_benefit);
                writer.WriteLine("$$$The Multiple_DAR_benefit is " + Multiple_DAR_benefit);


                for (int rr = 1; rr < r + 1; rr++)
                {
                    if (Multiple_DAR.GetValue(z_r[rr]) > 0.99)
                    {
                        par_z_r[rr] = 1;
                        
                        for (int j = W + 1; j < d1 + 1; j++)
                        {
                            if (a_rj[rr, j] == 1)
                            {
    
                                if (status_j[j] == 55)
                                {
                                    status_j[j] = 88;
                                    accepted_new_customer_roll[no_rolling] = accepted_new_customer_roll[no_rolling] + 1;
                                }
                            }
                        }
                    }
                    else
                    {
                        par_z_r[rr] = 0;
                    }
                }

                for (int j = W + 1; j < d1 + 1; j++)
                {
                    if (status_j[j] == 55)
                    {
                        status_j[j] = -88;
                         rejected_new_customer_roll[no_rolling] = rejected_new_customer_roll[no_rolling] + 1;
                    }
                }

                Console.WriteLine(" accepted_new_customer_roll[" + no_rolling + " ] = " + accepted_new_customer_roll[no_rolling]);
                Console.WriteLine(" rejected_new_customer_roll[" + no_rolling + " ] = " + rejected_new_customer_roll[no_rolling]);


                int NO_vehicle;
                int NO_vehicle_search = 0;
                int sub_route;
                int Next_vertex_to_find = 0;
                int searched_rr_for_zero = 0;

                for (int v = 1; v < W + 1; v++)
                {
                    first_cust_v[v] = 0;
                }

                NO_vehicle = 0;
                while (NO_vehicle_search < W)
                {
                    NO_vehicle_search = NO_vehicle_search + 1;
                    int ooo = -200;

                    sub_route = 0;
                    for (int rr = 1; rr < r + 1; rr++)
                    {
                        if (rr > searched_rr_for_zero && par_z_r[rr] == 1 && set_s[start_r[rr]] == -1)
                        {
                            ooo = 200;

                            NO_vehicle = vehicle_s[start_r[rr]];
                             searched_rr_for_zero = rr;
                            sub_route = sub_route + 1;
                            Next_vertex_to_find = n_location_s[end_r[rr]];

                            output_information_route_Discre0611_20210104_rr_rolling0202_ball_flexibledropoff_clearversion(rr);
                            output_information_route_Discre0611_20210104_rr_to_excel_flexibledropoff_rolling0202_0212_ball_matching_clearversion(NO_vehicle, sub_route, rr, next_resolve_end);

                            break;
                        }
                    }
                    while (Next_vertex_to_find != 2 * n + 1 && ooo == 200)
                    {
                        ooo = -200;
                        for (int rr = 1; rr < r + 1; rr++)
                        {
                            if (par_z_r[rr] == 1 && n_location_s[start_r[rr]] == Next_vertex_to_find)
                            {
                                ooo = 200;
                                sub_route = sub_route + 1;
                                Next_vertex_to_find = n_location_s[end_r[rr]];
                                if (Next_vertex_to_find != 2 * n + 1)
                                {
                                    output_information_route_Discre0611_20210104_rr_rolling0202_ball_flexibledropoff_clearversion(rr);
                                    output_information_route_Discre0611_20210104_rr_to_excel_flexibledropoff_rolling0202_0212_ball_matching_clearversion(NO_vehicle, sub_route, rr, next_resolve_end);
                                }
                                break;
                            }
                        }
                    }
                }
  
                Multiple_DAR.End();
            }
            else
            {
                Console.WriteLine("NO SOLUTION????");
                Console.Read();
                Multiple_DAR.End();
            }

            end_cplex0 = DateTime.Now;
            interval_cplex0 = end_cplex0.Subtract(begin_cplex0);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("CPLEX DONE! The Total Time is " + Convert.ToString(interval_cplex0.Hours) + ":" + Convert.ToString(interval_cplex0.Minutes) + ":" + Convert.ToString(interval_cplex0.Seconds));
            writer.WriteLine("*******************************************************************");
            writer.WriteLine("CPLEX DONE! The Total Time is " + Convert.ToString(interval_cplex0.Hours) + ":" + Convert.ToString(interval_cplex0.Minutes) + ":" + Convert.ToString(interval_cplex0.Seconds));
            Console.WriteLine("T_0 = " + T_0 + " customer = " + n + ";" + " Q = " + Q + "; apha = " + apha + ";");
            writer.WriteLine("T_0 = " + T_0 + " customer = " + n + ";" + " Q = " + Q + "; apha = " + apha + ";");

        }

      
   

  

        public const double rolling_unit = 10;
        public double[] EarliestTime_v = new double[W + 1];
        public int[] n_location_v = new int[W + 1];

        public int[] status_j = new int[n + 1];
        public int[] reposition_j = new int[n + 1];

        public const int P = 100000;

        public int[] No_release_subpath_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] time_start_subpath = new double[P];
        public double[] time_end_subpath = new double[P];

        public double[] detail_time_end_subpath = new double[P];
        public double[] detail_duration_subpath = new double[P];

        public double[] fare_subpath = new double[P];
        public double[] cost_subpath = new double[P];
        public double[] profit_subpath = new double[P];

        public int[] n_location_start_subpath = new int[P];
        public int[] n_location_end_subpath = new int[P];
        public int[,] a_subpathj = new int[P, n + 1];
        public int sub_path;

        public double[] accum_profit_v = new double[W + 1];
        public int[] vehicle_subpath = new int[P];

        public int[] stop_vehicle = new int[W + 1];

        public const int truen = 20000;

        public int[] status_truej = new int[truen + 1];

        public double[] aao_truej = new double[truen + 1];
        public double[] bbo_truej = new double[truen + 1];
        public double[] aad_truej = new double[truen + 1];
        public double[] bbd_truej = new double[truen + 1];

        public double[] deadline_truej = new double[truen + 1];

        public double[] q_truej = new double[2 * truen + 2];
        public double[] g_truej = new double[truen + 1];
        public double[] requesttime_truej = new double[truen + 1];

        public double[] release_aao_truej = new double[truen + 1];
        public double[] release_bbo_truej = new double[truen + 1];
        public double[] release_aad_truej = new double[truen + 1];
        public double[] release_bbd_truej = new double[truen + 1];


        public int[] release_aabbo_truej_data_candidate_no = new int[truen + 1];
        public int[] release_aabbd_truej_data_candidate_no = new int[truen + 1];

        public double[] release_start_time_truej = new double[2 * truen + 2];
        public double[] release_walk_time_truej = new double[2 * truen + 2];
        public double[] release_arrivehome_time_truej = new double[2 * truen + 2];
        public double[] release_arrivestop_time_truej = new double[2 * truen + 2];

        public int[] release_subpath_truej = new int[2 * truen + 2];

        public int[,] subpath_vehicle_stop = new int[W + 1, 2 * truen + 2];

        public int[,] truej_vehicle_stop = new int[W + 1, 2 * truen + 2];
        public int[] release_vehicle_truej = new int[2 * truen + 2];
        public int[] release_stop_seq_truej = new int[2 * truen + 2];

        public void Floating_target_results_output2(int no_vehicle, int all_control)
        {

            int[] FT_pickup_truej = new int[max_truej_untilnow + 1];
            int[] FT_dropoff_truej = new int[max_truej_untilnow + 1];
            int[] shareinpath_truej = new int[max_truej_untilnow + 1];

            string strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";


            if (all_control == 1)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_1.mdb";
            if (all_control == 2)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_2.mdb";
            if (all_control == 3)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_3.mdb";
            if (all_control == 4)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_4.mdb";
            if (all_control == 5)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_5.mdb";
            if (all_control == 6)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_6.mdb";
            if (all_control == 7)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_7.mdb";
            if (all_control == 8)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_8.mdb";
            if (all_control == 9)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_9.mdb";
            if (all_control == 10)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_10.mdb";
            if (all_control == 11)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_11.mdb";
            if (all_control == 12)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_12.mdb";
            if (all_control == 13)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_13.mdb";
            if (all_control == 14)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_14.mdb";
            if (all_control == 15)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_15.mdb";
            if (all_control == 16)
                strConnection += @"Data Source=|DataDirectory|\Output_Instance_file_v2_16.mdb";



            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str;
            int asqwqwq;

            str = "Delete from customer_coordinates";
            OleDbCommand myCommand1 = new OleDbCommand(str, objConnection);
            asqwqwq = myCommand1.ExecuteNonQuery();

            for (int v = 1; v < no_vehicle + 1; v++)
            {
                int truejjj = truej_vehicle_stop[v, 0];
                str = "INSERT INTO customer_coordinates (Vehicle, Stop, Customer, Passengers, Fare, Home_Intersection, Home_Lon, Home_Lat, FT_Intersection, FT_Lon, FT_Lat, Release_time, Walking_distance, Walking_time, Arrive_stop_time, Arrive_home_time, Deadline, sub_path, vehicle_sub_path, time_start, time_end, detail_time_end, detail_time_duration, sub_path_fare, sub_path_cost, sub_path_profit) VALUES (" + (v) + "," + (0) + "," + (-111) + "," + (0) + "," + (0) + ","
                + (truejjj - 100000) + "," + (aac_i[truejjj - 100000]) + "," + (bbc_i[truejjj - 100000]) + "," + (truejjj - 100000) + "," + (aac_i[truejjj - 100000]) + "," + (bbc_i[truejjj - 100000])
                + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111)
                + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + "," + (-111) + ")";
                OleDbCommand myCommand00 = new OleDbCommand(str, objConnection);
                asqwqwq = myCommand00.ExecuteNonQuery();

                for (int stop = 1; stop < stop_vehicle[v] + 1; stop++)
                {
                    truejjj = truej_vehicle_stop[v, stop];
                    if (truejjj > 100000)
                    {
                        int subp = subpath_vehicle_stop[v, stop];
                        int Home_Intersection = truejjj - 100000;
                        int FT_Intersection = truejjj - 100000;
                        str = "INSERT INTO customer_coordinates (Vehicle, Stop, Customer, Passengers, Fare, Home_Intersection, Home_Lon, Home_Lat, FT_Intersection, FT_Lon, FT_Lat, Release_time, Walking_distance, Walking_time, Arrive_stop_time, Arrive_home_time, Deadline, sub_path, vehicle_sub_path, time_start, time_end, detail_time_end, detail_time_duration, sub_path_fare, sub_path_cost, sub_path_profit) VALUES (" + (v) + "," + (stop) + "," + (truejjj) + "," + (0) + "," + (0) + ","
                        + Home_Intersection + "," + (aac_i[truejjj - 100000]) + "," + (bbc_i[truejjj - 100000]) + "," + FT_Intersection + "," + (aac_i[truejjj - 100000]) + "," + (bbc_i[truejjj - 100000])
                        + "," + (time_start_subpath[subp]) + "," + (-111) + "," + (-111) + "," + (detail_time_end_subpath[subp]) + "," + (-111) + "," + (-111)
                        + "," + (subp) + "," + (vehicle_subpath[subp]) + "," + (time_start_subpath[subp]) + "," + (time_end_subpath[subp]) + "," + (detail_time_end_subpath[subp]) + "," + (detail_duration_subpath[subp]) + "," + (fare_subpath[subp]) + "," + (cost_subpath[subp]) + "," + (profit_subpath[subp]) + ")";
                        OleDbCommand myCommand01 = new OleDbCommand(str, objConnection);
                        asqwqwq = myCommand01.ExecuteNonQuery();

                    }
                    else if (truejjj <= truen)
                    {
                        int subp = release_subpath_truej[truejjj];
                        int Home_Intersection = pickupi_match_j_no[truejjj, 1];
                        int FT_Intersection = pickupi_match_j_no[truejjj, release_aabbo_truej_data_candidate_no[truejjj]];

                        double WD = Real_Coordinate_Dis_Pedestrian_data_truej(truejjj, release_aabbo_truej_data_candidate_no[truejjj]);

                        if (WD > 0)
                        {
                            FT_pickup_truej[truejjj] = 1;
                        }
                        shareinpath_truej[truejjj] = shareno_subpath[subp];

                        str = "INSERT INTO customer_coordinates (Vehicle, Stop, Customer, Passengers, Fare, Home_Intersection, Home_Lon, Home_Lat, FT_Intersection, FT_Lon, FT_Lat, Release_time, Walking_distance, Walking_time, Arrive_stop_time, Arrive_home_time, Deadline, sub_path, vehicle_sub_path, time_start, time_end, detail_time_end, detail_time_duration, sub_path_fare, sub_path_cost, sub_path_profit) VALUES (" + (v) + "," + (stop) + "," + (truejjj) + "," + (q_truej[truejjj]) + "," + (g_truej[truejjj]) + ","
                             + Home_Intersection + "," + (aao_truej[truejjj]) + "," + (bbo_truej[truejjj]) + "," + FT_Intersection + "," + (release_aao_truej[truejjj]) + "," + (release_bbo_truej[truejjj])
                            + "," + (release_start_time_truej[truejjj]) + "," + (WD) + "," + (release_walk_time_truej[truejjj]) + "," + (release_arrivestop_time_truej[truejjj]) + "," + (-111) + "," + (-111) + "," + (release_subpath_truej[truejjj]) + "," + (vehicle_subpath[subp]) + "," + (time_start_subpath[subp]) + "," + (time_end_subpath[subp]) + "," + (detail_time_end_subpath[subp]) + "," + (detail_duration_subpath[subp]) + "," + (fare_subpath[subp]) + "," + (cost_subpath[subp]) + "," + (profit_subpath[subp]) + ")";
                        OleDbCommand myCommand01 = new OleDbCommand(str, objConnection);
                        asqwqwq = myCommand01.ExecuteNonQuery();
                    }
                    else
                    {
                        int subp = release_subpath_truej[truejjj];
                        int Home_Intersection = dropoffi_match_j_no[truejjj - truen, 1];
                        int FT_Intersection = dropoffi_match_j_no[truejjj - truen, release_aabbd_truej_data_candidate_no[truejjj - truen]];
                        double WD = Real_Coordinate_Dis_Pedestrian_data_truej(truejjj, release_aabbd_truej_data_candidate_no[truejjj - truen]);

                        if (WD > 0)
                        {
                            FT_dropoff_truej[truejjj - truen] = 1;
                        }

                        str = "INSERT INTO customer_coordinates (Vehicle, Stop, Customer, Passengers, Fare, Home_Intersection, Home_Lon, Home_Lat, FT_Intersection, FT_Lon, FT_Lat, Release_time, Walking_distance, Walking_time, Arrive_stop_time, Arrive_home_time, Deadline, sub_path, vehicle_sub_path, time_start, time_end, detail_time_end, detail_time_duration, sub_path_fare, sub_path_cost, sub_path_profit) VALUES (" + (v) + "," + (stop) + "," + (truejjj) + "," + (q_truej[truejjj]) + "," + (0) + ","
                             + Home_Intersection + "," + (aad_truej[truejjj - truen]) + "," + (bbd_truej[truejjj - truen]) + "," + FT_Intersection + "," + (release_aad_truej[truejjj - truen]) + "," + (release_bbd_truej[truejjj - truen])
                           + "," + (release_start_time_truej[truejjj]) + "," + (WD) + "," + (release_walk_time_truej[truejjj]) + "," + (release_arrivestop_time_truej[truejjj]) + "," + (release_arrivehome_time_truej[truejjj]) + "," + (deadline_truej[truejjj - truen]) + ","
                           + (release_subpath_truej[truejjj]) + "," + (vehicle_subpath[subp]) + "," + (time_start_subpath[subp]) + "," + (time_end_subpath[subp]) + "," + (detail_time_end_subpath[subp]) + "," + (detail_duration_subpath[subp]) + "," + (fare_subpath[subp]) + "," + (cost_subpath[subp]) + "," + (profit_subpath[subp]) + ")";
                        OleDbCommand myCommand01 = new OleDbCommand(str, objConnection);
                        asqwqwq = myCommand01.ExecuteNonQuery();
                    }
                }
            }

            str = "Delete from subpath_profits";
            OleDbCommand myCommand2 = new OleDbCommand(str, objConnection);
            asqwqwq = myCommand2.ExecuteNonQuery();

            for (int subp = 1; subp < sub_path + 1; subp++)
            {
                str = "INSERT INTO subpath_profits (sub_path, vehicle, time_start, time_end, detail_time_end, detail_time_duration, fare, cost, profit, shareno_subpath) VALUES (" + (subp) + "," + (vehicle_subpath[subp]) + "," + (time_start_subpath[subp]) + "," + (time_end_subpath[subp]) + "," + (detail_time_end_subpath[subp]) + "," + (detail_duration_subpath[subp]) + "," + (fare_subpath[subp]) + "," + (cost_subpath[subp]) + "," + (profit_subpath[subp]) + "," + shareno_subpath[subp] + ")";
                OleDbCommand myCommand02 = new OleDbCommand(str, objConnection);
                asqwqwq = myCommand02.ExecuteNonQuery();
            }

            str = "Delete from vehicle_profits";
            OleDbCommand myCommand4 = new OleDbCommand(str, objConnection);
            asqwqwq = myCommand4.ExecuteNonQuery();

            for (int v = 1; v < no_vehicle + 1; v++)
            {
                str = "INSERT INTO vehicle_profits (vehicle, accum_profit) VALUES (" + (v) + "," + (accum_profit_v[v]) + ")";
                OleDbCommand myCommand04 = new OleDbCommand(str, objConnection);
                asqwqwq = myCommand04.ExecuteNonQuery();
            }


            str = "Delete from customer_status";
            OleDbCommand myCommand3 = new OleDbCommand(str, objConnection);
            asqwqwq = myCommand3.ExecuteNonQuery();

            for (int truejjj = 1; truejjj < max_truej_untilnow + 1; truejjj++)
            {
                if (status_truej[truejjj] == 99)
                    str = "INSERT INTO customer_status (customer, status, sub_path, shareinpath_truej, FT_pickup_truej, pickup_WD, pickup_Home_Intersection, pickup_Home_Lon, pickup_Home_Lat, pickup_FT_Intersection, pickup_FT_Lon, pickup_FT_Lat, FT_dropoff_truej, dropoff_WD, dropoff_Home_Intersection, dropoff_Home_Lon, dropoff_Home_Lat, dropoff_FT_Intersection, dropoff_FT_Lon, dropoff_FT_Lat) VALUES (" + (truejjj) + "," + (1) + "," + release_subpath_truej[truejjj] + "," + shareinpath_truej[truejjj] + "," + FT_pickup_truej[truejjj] + "," + Real_Coordinate_Dis_Pedestrian_data_truej(truejjj, release_aabbo_truej_data_candidate_no[truejjj]) + "," + pickupi_match_j_no[truejjj, 1] + "," + (aao_truej[truejjj]) + "," + (bbo_truej[truejjj]) + "," + pickupi_match_j_no[truejjj, release_aabbo_truej_data_candidate_no[truejjj]] + "," + (release_aao_truej[truejjj]) + "," + (release_bbo_truej[truejjj]) + "," + FT_dropoff_truej[truejjj] + "," + Real_Coordinate_Dis_Pedestrian_data_truej(truejjj + truen, release_aabbd_truej_data_candidate_no[truejjj]) + "," + dropoffi_match_j_no[truejjj, 1] + "," + (aad_truej[truejjj]) + "," + (bbd_truej[truejjj]) + "," + dropoffi_match_j_no[truejjj, release_aabbd_truej_data_candidate_no[truejjj]] + "," + (release_aad_truej[truejjj]) + "," + (release_bbd_truej[truejjj]) + ")";
                else
                    str = "INSERT INTO customer_status (customer, status, sub_path, shareinpath_truej, FT_pickup_truej, pickup_WD, pickup_Home_Intersection, pickup_Home_Lon, pickup_Home_Lat, pickup_FT_Intersection, pickup_FT_Lon, pickup_FT_Lat, FT_dropoff_truej, dropoff_WD, dropoff_Home_Intersection, dropoff_Home_Lon, dropoff_Home_Lat, dropoff_FT_Intersection, dropoff_FT_Lon, dropoff_FT_Lat) VALUES (" + (truejjj) + "," + (0) + "," + release_subpath_truej[truejjj] + "," + shareinpath_truej[truejjj] + "," + FT_pickup_truej[truejjj] + "," + Real_Coordinate_Dis_Pedestrian_data_truej(truejjj, release_aabbo_truej_data_candidate_no[truejjj]) + "," + pickupi_match_j_no[truejjj, 1] + "," + (aao_truej[truejjj]) + "," + (bbo_truej[truejjj]) + "," + pickupi_match_j_no[truejjj, release_aabbo_truej_data_candidate_no[truejjj]] + "," + (release_aao_truej[truejjj]) + "," + (release_bbo_truej[truejjj]) + "," + FT_dropoff_truej[truejjj] + "," + Real_Coordinate_Dis_Pedestrian_data_truej(truejjj + truen, release_aabbd_truej_data_candidate_no[truejjj]) + "," + dropoffi_match_j_no[truejjj, 1] + "," + (aad_truej[truejjj]) + "," + (bbd_truej[truejjj]) + "," + dropoffi_match_j_no[truejjj, release_aabbd_truej_data_candidate_no[truejjj]] + "," + (release_aad_truej[truejjj]) + "," + (release_bbd_truej[truejjj]) + ")";
                OleDbCommand myCommand03 = new OleDbCommand(str, objConnection);
                asqwqwq = myCommand03.ExecuteNonQuery();
            }


            str = "Delete from rolling_iteration";
            OleDbCommand myCommand5 = new OleDbCommand(str, objConnection);
            asqwqwq = myCommand5.ExecuteNonQuery();

            for (int no_rolling = 1; no_rolling <= (int)Math.Ceiling(T_0 / rolling_unit); no_rolling++)
            {
                str = "INSERT INTO rolling_iteration (no_rolling, new_customer_roll, old_customer_roll, current_vehicle_roll, later_vehicle_roll, accepted_new_customer_roll, direct_rejected_new_customer_roll, rejected_new_customer_roll, release_customer_roll, release_vehicle_roll, idle_vehicle_roll, reposition_vehicle_roll, timespace_node_roll, timespace_route_roll, DP_CPU_time_roll) VALUES (" + (no_rolling) + "," + new_customer_roll[no_rolling] + "," + old_customer_roll[no_rolling] + "," + current_vehicle_roll[no_rolling] + "," + later_vehicle_roll[no_rolling] + "," + accepted_new_customer_roll[no_rolling] + "," + rejected_new_customer_roll[no_rolling] + "," + direct_rejected_new_customer_roll[no_rolling] + "," + release_customer_roll[no_rolling] + "," + release_vehicle_roll[no_rolling] + "," + idle_vehicle_roll[no_rolling] + "," + reposition_vehicle_roll[no_rolling] + "," + timespace_node_roll[no_rolling] + "," + timespace_route_roll[no_rolling] + "," + DP_CPU_time_roll[no_rolling] + ")";
                OleDbCommand myCommand05 = new OleDbCommand(str, objConnection);
                asqwqwq = myCommand05.ExecuteNonQuery();
            }


            objConnection.Close();
            Console.WriteLine("Output_data_finished!!!");
        }

       
        public int max_truej_untilnow;

        public double Dishere = 0;
        public double Accum_dishere = 0;
        public double Avg_dishere = 0;
        public double Max_dishere = 0;
        public int No_dishere = 0;
        public const int cand = 50;

        public int[] new_customer_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] old_customer_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] current_vehicle_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] later_vehicle_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] accepted_new_customer_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] rejected_new_customer_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] direct_rejected_new_customer_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] release_customer_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] release_vehicle_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] idle_vehicle_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] reposition_vehicle_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] timespace_node_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] timespace_route_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] DP_CPU_time_roll = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];

        public int[] shareno_subpath = new int[P];

        public void Rolling_Horizon_overall0202_0212_flexibledropoff_matching_clearversion_reposition_changej(int KKK, int no_vehicle, int all_control)
        {
            overall_control(all_control);

            double previous_resolve_start;
            double resolve_start;
            double resolve_end;
            double next_resolve_end;

            input_intersections();
            input_five_files_for_intersections_overall_data();

            input_vehicle_origin_at_intersections(no_vehicle);

            for (int v = 1; v < no_vehicle + 1; v++)
            {
                EarliestTime_v[v] = 0;
                n_location_v[v] = v + n;
            }
       
            sub_path = 0;
            No_release_subpath_roll[0] = 0;

            for (int initial_j = 0; initial_j < n + 1; initial_j++)
            {
                reposition_j[initial_j] = 0;
            }

            for (int no_rolling = 1; no_rolling <= (int)Math.Ceiling(T_0 / rolling_unit); no_rolling++)
            {
                previous_resolve_start = (no_rolling - 1) * rolling_unit;
                resolve_start = no_rolling * rolling_unit;
                resolve_end = (no_rolling + 1) * rolling_unit;
                next_resolve_end = (no_rolling + 2) * rolling_unit;
     
                Console.WriteLine("no_rolling = " + no_rolling + "; resolve_end = " + resolve_end);

                DateTime begin_total;
                DateTime end_total;
                TimeSpan interval_total;

                begin_total = DateTime.Now;

                input_ball_Discre0611_justinput_rolling_coordinates_matching_changej(previous_resolve_start, resolve_start, no_vehicle, Fle_pickup, Fle_dropoff);


                end_total = DateTime.Now;
                interval_total = end_total.Subtract(begin_total);
                Console.WriteLine("*******************************************************************");
                Console.WriteLine("DONE! The Total Time is " + interval_total.TotalSeconds);
                Console.WriteLine("input_ball_Discre0611_justinput_rolling_coordinates_matching_changej DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));


                for (int j = W + 1; j < d1 + 1; j++)
                {
                    if (status_j[j] == 55 || status_j[j] == 88)
                    {
                        Deadline_j[j] = Requesttime_j[j] + (1 + apha) * T_j[j];
                        if (status_j[j] == 55)
                            new_customer_roll[no_rolling] = new_customer_roll[no_rolling] + 1;
                        if (status_j[j] == 88)
                            old_customer_roll[no_rolling] = old_customer_roll[no_rolling] + 1;
                    }
                }
                Console.WriteLine(" new_customer_roll[" + no_rolling + " ] = " + new_customer_roll[no_rolling]);
                Console.WriteLine(" old_customer_roll[" + no_rolling + " ] = " + old_customer_roll[no_rolling]);


                for (int v = 1; v < no_vehicle + 1; v++)
                {
                    if (EarliestTime_v[v] <= resolve_end)
                    {
                        EarliestTime_v[v] = resolve_end;
                        current_vehicle_roll[no_rolling] = current_vehicle_roll[no_rolling] + 1;
                    }
                    else
                    {
                        later_vehicle_roll[no_rolling] = later_vehicle_roll[no_rolling] + 1;
                    }
                }
                Console.WriteLine(" current_vehicle_roll[" + no_rolling + " ] = " + current_vehicle_roll[no_rolling]);

              
                Multiple_DAR_whole_process_Restartversion_rolling_0202_ball_flexibledropoff_matching_clearversion(KKK, resolve_end, next_resolve_end, no_vehicle, no_rolling);

                No_release_subpath_roll[no_rolling] = sub_path;

                for (int pp = No_release_subpath_roll[no_rolling - 1] + 1; pp < No_release_subpath_roll[no_rolling] + 1; pp++)
                {
                    if (time_start_subpath[pp] < next_resolve_end)
                    {
                        accum_profit_v[vehicle_subpath[pp]] = accum_profit_v[vehicle_subpath[pp]] + profit_subpath[pp];

                        for (int j = 1; j < n + 1; j++)
                        {
                            if (a_subpathj[pp, j] == 1)
                            {
                                if (status_j[j] == 88)
                                {
                                    status_j[j] = 99;
                                    status_truej[true_j[j]] = 99;

                                    release_customer_roll[no_rolling] = release_customer_roll[no_rolling] + 1;

                                    shareno_subpath[pp] = shareno_subpath[pp] + 1;
                                    Dishere = Real_Coordinate_Dis_Car_data(n_location_v[vehicle_subpath[pp]], 1, j, 1);

                                    Accum_dishere = Accum_dishere + Dishere;
                                    No_dishere = No_dishere + 1;
                                    Avg_dishere = Accum_dishere / No_dishere;
                                    if (Dishere > Max_dishere)
                                        Max_dishere = Dishere;
                                }
                            }
                        }


                        if (EarliestTime_v[vehicle_subpath[pp]] < time_end_subpath[pp])
                        {
                            EarliestTime_v[vehicle_subpath[pp]] = time_end_subpath[pp];

                            if (n_location_v[vehicle_subpath[pp]] > n)
                            {
                                status_j[n_location_v[vehicle_subpath[pp]] - n] = 101;
                            }

                            n_location_v[vehicle_subpath[pp]] = n_location_end_subpath[pp];
                            status_j[n_location_v[vehicle_subpath[pp]] - n] = 100;
                            release_vehicle_roll[no_rolling] = release_vehicle_roll[no_rolling] + 1;
                        }
                    }
                }

                Console.WriteLine("release_customer_roll[" + no_rolling + "] = " + release_customer_roll[no_rolling]);
                Console.WriteLine("release_vehicle_roll[" + no_rolling + "] = " + release_vehicle_roll[no_rolling]);

                if (repositioning == 1)
                {
                    vehicle_repositioning2(resolve_end, no_rolling);
                }
            }
            for (int v = 1; v < no_vehicle + 1; v++)
            {
                Console.WriteLine(" >> Vehicle " + v + " earn profit = " + accum_profit_v[v]);
            }

            Floating_target_results_output2(no_vehicle, all_control);

        }

        public void vehicle_repositioning2(double resolve_end, int no_rolling) 
        {
            int idle_vehicle = 0;
            int reposition_vehicle = 0;
            int no_reposition_vehicle = 0;
            for (int v = 1; v < W + 1; v++)
            {
                if (EarliestTime_v[v] == resolve_end) 
                {
                    idle_vehicle = idle_vehicle + 1;

                    int compare_nno = 0;
                    int cal_no = 0;
                    for (int truejjj = max_truej_untilnow; truejjj > 0; truejjj--)
                    {
                        compare_nno = compare_nno + 1;
                        if (Real_Coordinate_Dis_Car_data_j_truej(n_location_v[v], 1, truejjj, 1) < Avg_dishere)
                        {
                            cal_no = cal_no + 1;
                        }
                        if (compare_nno == 1000)
                            break;
                    }

                    if (cal_no >= Math.Floor(repo_or_not_par * compare_nno))
                    {
                        no_reposition_vehicle = no_reposition_vehicle + 1;

                    }
                    else
                    {
                        reposition_vehicle = reposition_vehicle + 1;
     

                        int range = 60;

                    here_change_minute:;
                        int no_cand = 0;
                        int[] indexi_candi = new int[cand + 1];
                        int[] sumofcust_candi = new int[cand + 1];

                        while (no_cand == 0)
                        {
                    

                            for (int i = 1; i < Total_I + 1; i++)
                            {
                                double travel_time_here = Real_Coordinate_Time_Car_data_j_intersection(n_location_v[v], 1, i);

                                if (travel_time_here > range && travel_time_here <= range + 60) 
                                {
                                    no_cand = no_cand + 1;
                                    indexi_candi[no_cand] = i;
                                }

                                if (no_cand >= cand)
                                    break;
                            }
                            range = range + 60;
                        }
                        int max_indexi = -1;
                        int max_here = 0;

                        for (int iii = 1; iii < cand + 1; iii++)
                        {
                            int takei = indexi_candi[iii];
                            if (takei == 0)
                                break;

                            int compare_no = 0;
                            for (int truejjj = max_truej_untilnow; truejjj > 0; truejjj--)
                            {
                                compare_no = compare_no + 1;
                                if (Real_Coordinate_Dis_Car_data_intersection_truej(takei, truejjj, 1) < Avg_dishere)
                                {
                                    sumofcust_candi[iii] = sumofcust_candi[iii] + 1;
                                }
                                if (compare_no == 1000)
                                    break;
                            }
                            if (sumofcust_candi[iii] > max_here)
                            {
                                max_here = sumofcust_candi[iii];
                                max_indexi = indexi_candi[iii];
                            }
                        }

                        if (max_here == 0)
                        {
                            goto here_change_minute;
                        }
     
                        if (max_here > 0)
                        {
                          

                            sub_path = sub_path + 1;
                            vehicle_subpath[sub_path] = v;

                            time_start_subpath[sub_path] = resolve_end;
                            n_location_start_subpath[sub_path] = n_location_v[v];

                            detail_time_end_subpath[sub_path] = resolve_end + Real_Coordinate_Time_Car_data_j_intersection(n_location_v[v], 1, max_indexi);

                            detail_duration_subpath[sub_path] = Real_Coordinate_Time_Car_data_j_intersection(n_location_v[v], 1, max_indexi);

                            time_end_subpath[sub_path] = (int)Math.Ceiling((resolve_end + Real_Coordinate_Time_Car_data_j_intersection(n_location_v[v], 1, max_indexi)) / unit) * unit;


                            fare_subpath[sub_path] = 0;
                            cost_subpath[sub_path] = Real_Coordinate_Time_Car_data_j_intersection(n_location_v[v], 1, max_indexi) * cost;

                            profit_subpath[sub_path] = -Real_Coordinate_Time_Car_data_j_intersection(n_location_v[v], 1, max_indexi) * cost;

                            
                            stop_vehicle[v] = stop_vehicle[v] + 1;                           
                            subpath_vehicle_stop[v, stop_vehicle[v]] = sub_path;

                            int dddddd = n_location_v[v] - n;

                            true_j[dddddd] = 100000 + max_indexi;

                            aao_j[dddddd] = aac_i[max_indexi];
                            bbo_j[dddddd] = bbc_i[max_indexi];
                            aad_j[dddddd] = aac_i[max_indexi];
                            bbd_j[dddddd] = bbc_i[max_indexi];

                            n_location_end_subpath[sub_path] = dddddd + n;

                            truej_vehicle_stop[v, stop_vehicle[v]] = true_j[dddddd];

                            EarliestTime_v[v] = time_end_subpath[sub_path];

                            if (n_location_v[v] > n)
                            {
                                status_j[n_location_v[v] - n] = 101;
                            }

                            n_location_v[v] = n_location_end_subpath[sub_path];
                            status_j[n_location_v[v] - n] = 100;

                            reposition_j[n_location_v[v] - n] = -1;
                        }
                    }
                }
            }
            Console.WriteLine("repo_or_not_par = " + repo_or_not_par + " : " + reposition_vehicle + " out of " + idle_vehicle + " DO reposition <=> " + no_reposition_vehicle + " out of " + idle_vehicle + " DONOT reposition");
            idle_vehicle_roll[no_rolling] = idle_vehicle;
            reposition_vehicle_roll[no_rolling] = reposition_vehicle;
        }
     
      
        public void Multiple_DAR_whole_process_Restartversion_rolling_0202_ball_flexibledropoff_matching_clearversion(int KKK, double resolve_end, double next_resolve_end, int vehicle_no, int no_rolling)
        {
            r = 0; 
            r_dummy = 0; 
            s = 0;

            DateTime begin_total0;
            DateTime end_total0;
            TimeSpan interval_total0;

            begin_total0 = DateTime.Now;

            Console.WriteLine("Multiple_DAR_whole_process_Restartversion_rolling_0202_ball_flexibledropoff_matching_clearversion >> T_0 = " + T_0 + " customer = " + n + ";" + " Q = " + Q + "; apha = " + apha + ";");
        

            create_time_space_within_timewindow_for_rolling_0202_ball_flexibledropoff_matching_changej_0405(KKK, resolve_end, no_rolling);
            
            Console.WriteLine("TOTAL  TOTAL  TOTAL  TOTAL  TOTAL  r = " + r);

            timespace_route_roll[no_rolling] = r;
            Console.WriteLine("timespace_route_roll[" + no_rolling + "] = " + timespace_route_roll[no_rolling]);


            DateTime begin_cplex;
            DateTime end_cplex;
            TimeSpan interval_cplex;

            begin_cplex = DateTime.Now;

            deterministic_fomulation_Discre_rolling_0202_ball_flexibledropoff_matching_clearversion_changej_improve(vehicle_no, next_resolve_end, no_rolling);

            end_cplex = DateTime.Now;
            interval_cplex = end_cplex.Subtract(begin_cplex);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("ALL CPLEX DONE! The Total Time is " + Convert.ToString(interval_cplex.Hours) + ":" + Convert.ToString(interval_cplex.Minutes) + ":" + Convert.ToString(interval_cplex.Seconds));
            writer.WriteLine("*******************************************************************");
            writer.WriteLine("ALL CPLEX DONE! The Total Time is " + Convert.ToString(interval_cplex.Hours) + ":" + Convert.ToString(interval_cplex.Minutes) + ":" + Convert.ToString(interval_cplex.Seconds));
            Console.WriteLine("T_0 = " + T_0 + " customer = " + n + ";" + " Q = " + Q + "; apha = " + apha + ";");
            writer.WriteLine("T_0 = " + T_0 + " customer = " + n + ";" + " Q = " + Q + "; apha = " + apha + ";");


            end_total0 = DateTime.Now;
            interval_total0 = end_total0.Subtract(begin_total0);
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("DONE! The Total Time is " + Convert.ToString(interval_total0.Hours) + ":" + Convert.ToString(interval_total0.Minutes) + ":" + Convert.ToString(interval_total0.Seconds));
            Console.WriteLine("DONE! The Total Time is " + Convert.ToInt32(interval_total0.TotalSeconds));
            writer.WriteLine("*******************************************************************");
            writer.WriteLine("DONE! The Total Time is " + Convert.ToString(interval_total0.Hours) + ":" + Convert.ToString(interval_total0.Minutes) + ":" + Convert.ToString(interval_total0.Seconds));
            writer.WriteLine("DONE! The Total Time is " + Convert.ToInt32(interval_total0.TotalSeconds));
        }
   
      
        public void single_target_optimization_ball_Discre0611_flexibledropoff_changej(int iii, double Start_time)
        {
            
            Three_points_updating_flexibledropoff_changej(iii);
            Case_for_Ball_dis_Discre0611_flexibledropoff(iii, Start_time);

            if (j_3 <= n)
            {
                aa_ln[iii] = aav_i[iii];
                bb_ln[iii] = bbv_i[iii];

                aabb_ln_data[iii] = aavbbv_i_data[iii];
                aabb_ln_data_candidate_no[iii] = aavbbv_i_data_candidate_no[iii];
              
                if (Real_Coordinate_Time_Pedestrian_data(aavbbv_i_data[iii], aavbbv_i_data_candidate_no[iii]) > T_floating[iii - 1] + Real_Coordinate_Time_Car_data(xy_1_data, xy_1_data_candidate_no, aavbbv_i_data[iii], aavbbv_i_data_candidate_no[iii]))
                {                   
                    T_floating[iii] = Real_Coordinate_Time_Pedestrian_data(aavbbv_i_data[iii], aavbbv_i_data_candidate_no[iii]);
                }
                else
                {                  
                    T_floating[iii] = T_floating[iii - 1] + Real_Coordinate_Time_Car_data(aavbbv_i_data[iii - 1], aavbbv_i_data_candidate_no[iii - 1], aavbbv_i_data[iii], aavbbv_i_data_candidate_no[iii]);
                }

            }

            if (j_3 > n && feasibility == 1)
            {
                aa_ln[iii] = aav_i[iii];
                bb_ln[iii] = bbv_i[iii];

                aabb_ln_data[iii] = aavbbv_i_data[iii];
                aabb_ln_data_candidate_no[iii] = aavbbv_i_data_candidate_no[iii];
                
                T_floating[iii] = T_floating[iii - 1] + Real_Coordinate_Time_Car_data(aavbbv_i_data[iii - 1], aavbbv_i_data_candidate_no[iii - 1], aavbbv_i_data[iii], aavbbv_i_data_candidate_no[iii]);
            }

            if (j_3 > n && feasibility == -1)
            {
                aa_ln[iii] = aav_i[iii];
                bb_ln[iii] = bbv_i[iii];

                aabb_ln_data[iii] = aavbbv_i_data[iii];
                aabb_ln_data_candidate_no[iii] = aavbbv_i_data_candidate_no[iii];
                T_floating[iii] = -1000;
            }
        }

    
        int feasibility = 0;
        public void Case_for_Ball_dis_Discre0611_flexibledropoff(int iii, double Start_time)
        {
            double check_x, check_y;

            double check_t = 0;
            double min_here = int.MaxValue;
            int best_c = -1;

            feasibility = -1;

            for (int c = 1; c < discre_no + 1; c++)
            {
                check_x = candidate_aam_jc[j_3, c];
                check_y = candidate_bbm_jc[j_3, c];

              
                if (Real_Coordinate_Dis_Pedestrian_data(j_3, c) <= lw_j[j_3] + 0.0001)
                {
                    if (j_3 <= n)
                    {
                        feasibility = 1;
                        if (Real_Coordinate_Time_Pedestrian_data(j_3, c) > T_arrivex1 + Real_Coordinate_Time_Car_data(xy_1_data, xy_1_data_candidate_no, j_3, c))
                        {
                            check_t = Real_Coordinate_Time_Pedestrian_data(j_3, c) + Real_Coordinate_Time_Car_data(j_3, c, xy_2_data, xy_2_data_candidate_no);
                        }
                        else
                        {
                            check_t = T_arrivex1 + Real_Coordinate_Time_Car_data(xy_1_data, xy_1_data_candidate_no, j_3, c) + Real_Coordinate_Time_Car_data(j_3, c, xy_2_data, xy_2_data_candidate_no);
                        }
                    }
                    if (j_3 > n)
                    {
                        if (T_arrivex1 + Real_Coordinate_Time_Car_data(xy_1_data, xy_1_data_candidate_no, j_3, c) + Real_Coordinate_Time_Pedestrian_data(j_3, c) <= Deadline_j[j_3 - n] - Start_time)
                        {
                            check_t = T_arrivex1 + Real_Coordinate_Time_Car_data(xy_1_data, xy_1_data_candidate_no, j_3, c);
                            feasibility = 1;
                        }
                        else
                        {
                            check_t = 9999999999;
                        }
                    }


                    if (check_t < min_here)
                    {
                        min_here = check_t;
                        best_c = c;


                    }
                }
            }

            if (feasibility == 1)
            {
                case_ln[iii] = 1000 + best_c;
                aavbbv_i_data[iii] = j_3;
                aavbbv_i_data_candidate_no[iii] = best_c;

            }
            if (feasibility == -1)
            {
                case_ln[iii] = -1000;
                aavbbv_i_data[iii] = -1000;
                aavbbv_i_data_candidate_no[iii] = -1000;
            }

        }

    
        void Bit()
        {
            int ooo = 0;
            int ipow;
            ipow = 0;
            for (int j = 1; j < n + 1; j++)
            {
                lwrd_j[j] = (int)((j - 1) / 30) + 1;
                lbit_j[j] = (int)Math.Pow(2, ipow);
                ooo = ooo + lbit_j[j];
                ipow = ipow + 1;
                if (ipow == 30)
                {
                    ipow = 0;
                    ooo = 0;
                }
            }

            nwords_2 = (int)(d1 / 30) + 1;         
        }

       
        void Dominance_rule_nwords1(int ll)
        {

            if (d_l[ll] == 1)
            {
                for (int l = 0; l < Lall[0]; l++)
                {
                    if (d_l[l] == 1 && l != ll)
                    {
                        if (n_l[l] == n_l[ll])
                        {
                            if (T_l[l] <= T_l[ll] && pipipi_l[l] <= pipipi_l[ll])
                            {

                                int kkk = 0;
                                for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                {
                                    if ((preset1_ll[l, iw] | preset1_ll[ll, iw]) == preset1_ll[ll, iw])
                                    {
                                        kkk = kkk + 1;
                                        if ((preset2_ll[l, iw] | preset2_ll[ll, iw]) == preset2_ll[ll, iw])
                                        {
                                            kkk = kkk + 1;
                                        }
                                        else
                                        {
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                                if (kkk == (nwords_2 + 1 - nwords_1) * 2)
                                {
                                    d_l[ll] = 0;
                                    break;
                                }
                              
                            }

                            if (T_l[ll] <= T_l[l] && pipipi_l[ll] <= pipipi_l[l])
                            {
                               
                                int ggg = 0;
                                for (int iw = nwords_1; iw < nwords_2 + 1; iw++)
                                {
                                    if ((preset1_ll[ll, iw] | preset1_ll[l, iw]) == preset1_ll[l, iw])
                                    {
                                        ggg = ggg + 1;
                                        if ((preset2_ll[ll, iw] | preset2_ll[l, iw]) == preset2_ll[l, iw])
                                        {
                                            ggg = ggg + 1;
                                        }
                                        else
                                        {
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                                if (ggg == (nwords_2 + 1 - nwords_1) * 2)
                                {
                                    d_l[l] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
       

        public int[,] x_ij = new int[2 * nprime + 2, 2 * n + 2];

        public double x_1, y_1, T_arrivex1, x_2, y_2, x_3, y_3, r_3, s_3;
        public int j_3;

        public double[] aav_i = new double[2 * nprime + 2];
        public double[] bbv_i = new double[2 * nprime + 2];

        public double[] previous_aav_i = new double[2 * nprime + 2];
        public double[] previous_bbv_i = new double[2 * nprime + 2];

        double[] T_floating = new double[2 * nprime + 2];
        public double[] T_updating = new double[L];

        public double[] T_updatingR = new double[R + 1];

        public double x_0, y_0;



        public int xy_0_data, xy_1_data, xy_2_data, xy_3_data;
        public int xy_0_data_candidate_no, xy_1_data_candidate_no, xy_2_data_candidate_no, xy_3_data_candidate_no;

        void Three_points_updating_flexibledropoff_changej(int iii)
        {
            if (iii - 2 >= previous_drop_off_location)
            {
                x_0 = aav_i[iii - 2];
                y_0 = bbv_i[iii - 2];

                xy_0_data = aavbbv_i_data[iii - 2];
                xy_0_data_candidate_no = aavbbv_i_data_candidate_no[iii - 2];
            }

            x_1 = aav_i[iii - 1];
            y_1 = bbv_i[iii - 1];

            xy_1_data = aavbbv_i_data[iii - 1];
            xy_1_data_candidate_no = aavbbv_i_data_candidate_no[iii - 1];

            T_arrivex1 = T_floating[iii - 1];

            x_2 = aav_i[iii + 1];
            y_2 = bbv_i[iii + 1];

            xy_2_data = aavbbv_i_data[iii + 1];
            xy_2_data_candidate_no = aavbbv_i_data_candidate_no[iii + 1];

            for (int jjj = W + 1; jjj < d1 + 1; jjj++) 
            {
                if (x_ij[iii, jjj] == 1)
                {
                    j_3 = jjj;
                    x_3 = aao_j[jjj];
                    y_3 = bbo_j[jjj];

                    xy_3_data = jjj;
                    xy_3_data_candidate_no = 1;

                    r_3 = lw_j[jjj];
                    s_3 = sw_j[jjj];
                    break;
                }
            }

            for (int jjj = n + 1; jjj < (n + d1) + 1; jjj++) 
            {
                if (x_ij[iii, jjj] == 1)
                {
                    j_3 = jjj;
                    x_3 = aad_j[jjj - n];
                    y_3 = bbd_j[jjj - n];

                    xy_3_data = jjj;
                    xy_3_data_candidate_no = 1;

                    r_3 = lw_j[jjj - n];
                    s_3 = sw_j[jjj - n];
                    break;
                }
            }

        }

        public double mm_eps = 0.1;
        public int pppppp = 0;

        public int[] aavbbv_i_data = new int[2 * nprime + 2];
        public int[] aavbbv_i_data_candidate_no = new int[2 * nprime + 2];

        public int[] previous_aavbbv_i_data = new int[2 * nprime + 2];
        public int[] previous_aavbbv_i_data_candidate_no = new int[2 * nprime + 2];

        public int[] aabb_ln_data = new int[100];
        public int[] aabb_ln_data_candidate_no = new int[100];
        void One_route_updating_combined_ball_Discre0611_flexibledropoff_changej(int ll, double start_time)
        {
            aav_i[0] = aao_j[0];
            bbv_i[0] = bbo_j[0];

            aavbbv_i_data[0] = 0;
            aavbbv_i_data_candidate_no[0] = 1;

            for (int xx = 0; xx < 2 * nprime + 2; xx++)
            {
                for (int yy = 0; yy < 2 * n + 2; yy++)
                {
                    x_ij[xx, yy] = 0;
                }
            }

            x_ij[S_l[ll], n_l[ll]] = 1;
           
            aav_i[S_l[ll]] = aad_j[n_l[ll] - n];
            bbv_i[S_l[ll]] = bbd_j[n_l[ll] - n];

            aavbbv_i_data[S_l[ll]] = n_l[ll];
            aavbbv_i_data_candidate_no[S_l[ll]] = 1;

            int pre_label = ll;
            int pre_vertex = 0;
            for (int ii = 0; ii < S_l[ll]; ii++)
            {
                pre_label = prel_l[pre_label];
                pre_vertex = n_l[pre_label];

                previous_drop_off_location = 0;
                previous_drop_off_point = 0;
                T_arriving_previous_drop_off = 0;
                g_arriving_previous_drop_off = 0;

                if (pre_vertex <= n && pre_vertex > 0)
                {
                    x_ij[S_l[ll] - 1 - ii, pre_vertex] = 1;
                    aav_i[S_l[ll] - 1 - ii] = aao_j[pre_vertex];
                    bbv_i[S_l[ll] - 1 - ii] = bbo_j[pre_vertex];

                    aavbbv_i_data[S_l[ll] - 1 - ii] = pre_vertex;
                    aavbbv_i_data_candidate_no[S_l[ll] - 1 - ii] = 1;
                }
                else if (pre_vertex == 0) 
                {
                    x_ij[S_l[ll] - 1 - ii, pre_vertex] = 1;
                    aav_i[S_l[ll] - 1 - ii] = aao_j[pre_vertex];
                    bbv_i[S_l[ll] - 1 - ii] = bbo_j[pre_vertex];

                    aavbbv_i_data[S_l[ll] - 1 - ii] = pre_vertex;
                    aavbbv_i_data_candidate_no[S_l[ll] - 1 - ii] = 1;

                    previous_drop_off_point = pre_vertex;
                    previous_drop_off_location = S_l[pre_label];
                    T_arriving_previous_drop_off = T_l[pre_label];
                    g_arriving_previous_drop_off = -(pipipi_l[pre_label] - T_l[pre_label] * cost);
                }
                else
                {
                    x_ij[S_l[ll] - 1 - ii, pre_vertex] = 1;
                    previous_drop_off_point = pre_vertex;
                    previous_drop_off_location = S_l[pre_label];
                    T_arriving_previous_drop_off = T_l[pre_label];
                    g_arriving_previous_drop_off = -(pipipi_l[pre_label] - T_l[pre_label] * cost);

                    aav_i[S_l[ll] - 1 - ii] = aad_j[pre_vertex - n];
                    bbv_i[S_l[ll] - 1 - ii] = bbd_j[pre_vertex - n];

                    aavbbv_i_data[S_l[ll] - 1 - ii] = pre_vertex;
                    aavbbv_i_data_candidate_no[S_l[ll] - 1 - ii] = 1;
                }
            }

            double T_previous = int.MaxValue;
            T_floating[previous_drop_off_location] = T_arriving_previous_drop_off;

            int iteration = 0;
            T_updating[ll] = int.MinValue;
            while (T_updating[ll] + eps < T_previous)
            {
                if (iteration > 0) T_previous = T_updating[ll];
                iteration = iteration + 1;

                for (int iii = previous_drop_off_location + 1; iii < S_l[ll]; iii++)
                {
                    single_target_optimization_ball_Discre0611_flexibledropoff_changej(iii, start_time);
                    if (feasibility == -1)
                    {
                        goto End;
                    }
                }

                T_updating[ll] = T_floating[S_l[ll] - 1] + Real_Coordinate_Time_Car_data(aavbbv_i_data[S_l[ll] - 1], aavbbv_i_data_candidate_no[S_l[ll] - 1], n_l[ll], 1);
   
                aa_ln[S_l[ll]] = aad_j[n_l[ll] - n];
                bb_ln[S_l[ll]] = bbd_j[n_l[ll] - n];
                case_ln[S_l[ll]] = 8;

                aabb_ln_data[S_l[ll]] = n_l[ll];
                aabb_ln_data_candidate_no[S_l[ll]] = 1;
            }


            int k = S_l[ll] - 1;

            int no = 0;

            while (k > previous_drop_off_location)
            {
                no = no + 1;

                if (no > 100)
                    break;

                previous_aav_i[k] = aav_i[k];
                previous_bbv_i[k] = bbv_i[k];

                previous_aavbbv_i_data[k] = aavbbv_i_data[k];
                previous_aavbbv_i_data_candidate_no[k] = aavbbv_i_data_candidate_no[k];

                aav_i[k] = 100;
                bbv_i[k] = 100;
                
                single_target_optimization_ball_Discre0611_flexibledropoff_changej(k, start_time);
                if (feasibility == -1)
                {
                    Console.WriteLine("Position 2, goto end; something wrong");
                    goto End;
                }

                if (Real_Coordinate_Dis_Car_data(previous_aavbbv_i_data[k], previous_aavbbv_i_data_candidate_no[k], aavbbv_i_data[k], aavbbv_i_data_candidate_no[k]) <= mm_eps)
                {         
                    k = k - 1;
                }
                else
                {
                    if (k == S_l[ll] - 1)
                    {
                        k = k - 1;               
                    }
                    else
                    {
                        for (int kk = k + 1; kk <= S_l[ll] - 1; kk++)
                        {
                            single_target_optimization_ball_Discre0611_flexibledropoff_changej(kk, start_time);
                            if (feasibility == -1)
                            {
                                Console.WriteLine("Position 3, goto end; something wrong");
                                goto End;
                            }
                        }
                        k = S_l[ll] - 1;
                    }
                }
            }

           
            T_updating[ll] = T_floating[S_l[ll] - 1] + Real_Coordinate_Time_Car_data(aavbbv_i_data[S_l[ll] - 1], aavbbv_i_data_candidate_no[S_l[ll] - 1], n_l[ll], 1);

            aa_ln[S_l[ll]] = aad_j[n_l[ll] - n];
            bb_ln[S_l[ll]] = bbd_j[n_l[ll] - n];

            aabb_ln_data[S_l[ll]] = n_l[ll];
            aabb_ln_data_candidate_no[S_l[ll]] = 1;

            case_ln[S_l[ll]] = 8;


            T_l[ll] = T_updating[ll];

            double g = 0;
            for (int i = previous_drop_off_location + 1; i < S_l[ll] + 1; i++)
            {
                for (int j = W + 1; j < d1 + 1; j++)
                {
                    g = g + g_j[j] * x_ij[i, j];
                }
            }

            g = g + g_arriving_previous_drop_off;

            pipipi_l[ll] = -(g - T_l[ll] * cost);

        End:;
        }

        void One_route_updating_combined_ball_Discre0611_20210104_rr_flexibledropoff(int rr, int last_pre_vertex)
        {
            aav_i[0] = aao_j[0];
            bbv_i[0] = bbo_j[0];

            aavbbv_i_data[0] = 0;
            aavbbv_i_data_candidate_no[0] = 1;

            for (int xx = 0; xx < 2 * nprime + 2; xx++)
            {
                for (int yy = 0; yy < 2 * n + 2; yy++)
                {
                    x_ij[xx, yy] = 0;
                }
            }

            x_ij[seq_rj[rr, last_pre_vertex] - 999, last_pre_vertex] = 1;
            aav_i[seq_rj[rr, last_pre_vertex] - 999] = aad_j[last_pre_vertex - n];
            bbv_i[seq_rj[rr, last_pre_vertex] - 999] = bbd_j[last_pre_vertex - n];

            aavbbv_i_data[seq_rj[rr, last_pre_vertex] - 999] = last_pre_vertex;
            aavbbv_i_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999] = 1;

            int pre_vertex = -1;
            for (int ii = 0; ii < seq_rj[rr, last_pre_vertex] - 999; ii++)
            {
                for (int jjjj = 0; jjjj < 2 * n + 1; jjjj++)
                {
                    if (seq_rj[rr, jjjj] == seq_rj[rr, last_pre_vertex] - ii - 1)
                    {
                        pre_vertex = jjjj;
                        break;
                    }
                }

                previous_drop_off_location = 0;
                previous_drop_off_point = 0;
                T_arriving_previous_drop_off = 0;
                g_arriving_previous_drop_off = 0;

                if (pre_vertex <= n && pre_vertex > 0)
                {
                    x_ij[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii, pre_vertex] = 1;
                   
                    aav_i[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = aao_j[pre_vertex];
                    bbv_i[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = bbo_j[pre_vertex];

                    aavbbv_i_data[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = pre_vertex;
                    aavbbv_i_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = 1;

                }
                else if (pre_vertex == 0)
                {
                    x_ij[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii, pre_vertex] = 1;
                    aav_i[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = aao_j[pre_vertex];
                    bbv_i[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = bbo_j[pre_vertex];

                    aavbbv_i_data[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = pre_vertex;
                    aavbbv_i_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = 1;

                    previous_drop_off_point = pre_vertex;
                    previous_drop_off_location = seq_rj[rr, pre_vertex] - 999;
                    T_arriving_previous_drop_off = T_rj[rr, pre_vertex];
                    g_arriving_previous_drop_off = -(pipipi_rj[rr, pre_vertex] - T_rj[rr, pre_vertex] * cost);
                }
                else
                {
                    x_ij[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii, pre_vertex] = 1;
                    previous_drop_off_point = pre_vertex;
                    previous_drop_off_location = seq_rj[rr, pre_vertex] - 999;
                    T_arriving_previous_drop_off = T_rj[rr, pre_vertex];
                    g_arriving_previous_drop_off = -(pipipi_rj[rr, pre_vertex] - T_rj[rr, pre_vertex] * cost);

                    aav_i[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = aad_j[pre_vertex - n];
                    bbv_i[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = bbd_j[pre_vertex - n];

                    aavbbv_i_data[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = pre_vertex;
                    aavbbv_i_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999 - 1 - ii] = 1;
                }
            }
            double T_previous = int.MaxValue;

            T_floating[previous_drop_off_location] = T_arriving_previous_drop_off;


            int iteration = 0;
            T_updatingR[rr] = int.MinValue;
            while (T_updatingR[rr] + eps < T_previous)
            {
                if (iteration > 0) T_previous = T_updatingR[rr];
                iteration = iteration + 1;

                for (int iii = previous_drop_off_location + 1; iii < seq_rj[rr, last_pre_vertex] - 999; iii++)
                {
                    single_target_optimization_ball_Discre0611_flexibledropoff_changej(iii, time_s[start_r[rr]] * unit);
                }

                T_updatingR[rr] = T_floating[seq_rj[rr, last_pre_vertex] - 999 - 1] + Real_Coordinate_Time_Car_data(aavbbv_i_data[seq_rj[rr, last_pre_vertex] - 999 - 1], aavbbv_i_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999 - 1], last_pre_vertex, 1);
                aa_ln[seq_rj[rr, last_pre_vertex] - 999] = aad_j[last_pre_vertex - n];
                bb_ln[seq_rj[rr, last_pre_vertex] - 999] = bbd_j[last_pre_vertex - n];
                aabb_ln_data[seq_rj[rr, last_pre_vertex] - 999] = last_pre_vertex;
                aabb_ln_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999] = 1;
                case_ln[seq_rj[rr, last_pre_vertex] - 999] = 8;
            }

            int k = seq_rj[rr, last_pre_vertex] - 999 - 1;
            int no = 0;

            while (k > previous_drop_off_location)
            {
                no = no + 1;

                if (no > 100)
                    break;

                previous_aav_i[k] = aav_i[k];
                previous_bbv_i[k] = bbv_i[k];

                previous_aavbbv_i_data[k] = aavbbv_i_data[k];
                previous_aavbbv_i_data_candidate_no[k] = aavbbv_i_data_candidate_no[k];

                aav_i[k] = 100;
                bbv_i[k] = 100;
              
                single_target_optimization_ball_Discre0611_flexibledropoff_changej(k, time_s[start_r[rr]] * unit);
        
                if (Real_Coordinate_Dis_Car_data(previous_aavbbv_i_data[k], previous_aavbbv_i_data_candidate_no[k], aavbbv_i_data[k], aavbbv_i_data_candidate_no[k]) <= mm_eps)
                {                  
                    k = k - 1;
                }
                else
                {
                    if (k == seq_rj[rr, last_pre_vertex] - 999 - 1)
                    {
                        k = k - 1;             
                    }
                    else
                    {
                        for (int kk = k + 1; kk <= seq_rj[rr, last_pre_vertex] - 999 - 1; kk++)
                        {
                            single_target_optimization_ball_Discre0611_flexibledropoff_changej(kk, time_s[start_r[rr]] * unit);
                        }
                        k = seq_rj[rr, last_pre_vertex] - 999 - 1;   
                    }
                }
            }

            T_updatingR[rr] = T_floating[seq_rj[rr, last_pre_vertex] - 999 - 1] + Real_Coordinate_Time_Car_data(aavbbv_i_data[seq_rj[rr, last_pre_vertex] - 999 - 1], aavbbv_i_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999 - 1], last_pre_vertex, 1);
            aa_ln[seq_rj[rr, last_pre_vertex] - 999] = aad_j[last_pre_vertex - n];
            bb_ln[seq_rj[rr, last_pre_vertex] - 999] = bbd_j[last_pre_vertex - n];
            aabb_ln_data[seq_rj[rr, last_pre_vertex] - 999] = last_pre_vertex;
            aabb_ln_data_candidate_no[seq_rj[rr, last_pre_vertex] - 999] = 1;
            case_ln[seq_rj[rr, last_pre_vertex] - 999] = 8;

            T_rj[rr, last_pre_vertex] = (float)T_updatingR[rr];
            double g = 0;
            for (int i = previous_drop_off_location + 1; i < seq_rj[rr, last_pre_vertex] - 999 + 1; i++)
            {
                for (int j = 1; j < n + 1; j++)
                {
                    g = g + g_j[j] * x_ij[i, j];
                }
            }
            g = g + g_arriving_previous_drop_off;
            pipipi_rj[rr, last_pre_vertex] = (float)(-(g - T_rj[rr, last_pre_vertex] * cost));
        }


        //===========================================================================================================================


        public double[] Profit_control = new double[20];
        public int[] Total_customer_control = new int[20];
        public int[] Acceptance_control = new int[20];
        public double[] Acceptance_rate_control = new double[20];
        public int[] All_Floating_control = new int[20];
        public int[] FT_Pickup_control = new int[20];
        public int[] FT_Dropoff_control = new int[20];
        public double[] All_Floating_rate_control = new double[20];
        public double[] FT_Pickup_rate_control = new double[20];
        public double[] FT_Dropoff_rate_control = new double[20];
        public int[] Total_unreposition_path_control = new int[20];
        public int[] Shared_path_control = new int[20];
        public double[] Shared_path_rate_control = new double[20];

        public double[] Total_fare_control = new double[20];
        public double[] Total_cost_customer_control = new double[20];
        public double[] Total_cost_reposition_control = new double[20];
        public double[] Total_profit_control = new double[20];

        public double[] Duration_not_shared_subpath_control = new double[20];
        public double[] Duration_shared_subpath_control = new double[20];
        public double[] Duration_all_subpath_control = new double[20];


        public double[] new_customer_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] old_customer_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] current_vehicle_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] later_vehicle_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] accepted_new_customer_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] rejected_new_customer_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] direct_rejected_new_customer_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] release_customer_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] release_vehicle_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] idle_vehicle_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public double[] reposition_vehicle_avg_control = new double[(int)Math.Ceiling(T_0 / rolling_unit) + 1];

        public int[] timespace_node_1_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] timespace_node_2_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] timespace_node_3_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];

        public int[] timespace_route_1_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] timespace_route_2_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] timespace_route_3_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];

        public int[] DP_CPU_time_1_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] DP_CPU_time_2_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];
        public int[] DP_CPU_time_3_control = new int[(int)Math.Ceiling(T_0 / rolling_unit) + 1];

     
    
        public void Floating_target_results_Excel_Table0921_ALL_220421(int all_control, int TYPE)
        {
            int date = TYPE - 1200;
            input_all_fare_0921(date);

            Floating_target_results_Excel_Table0921_220421(all_control, TYPE);

            Console.Write(pickup_WD_total_control[all_control] + "    ");
            Console.WriteLine();
            Console.Write(dropoff_WD_total_control[all_control] + "    ");
            Console.WriteLine();
            Console.Write((pickup_WD_total_control[all_control] + dropoff_WD_total_control[all_control]) + "    ");
            Console.WriteLine();
            Console.Write(Fare_FT_control[all_control] + "    ");
            Console.WriteLine();


            input_all_Q_0921_220422(date);
            Floating_target_results_Excel_Table0921_220421_220422(all_control, TYPE);
            Console.Write(accept_Q_times_request_pickup_FT_total_control[all_control] + "  = accept_Q_times_request_pickup_FT_total_control  ");
            Console.WriteLine();
            Console.Write(accept_Q_times_request_dropoff_FT_total_control[all_control] + "   = accept_Q_times_request_dropoff_FT_total_control ");
            Console.WriteLine();
            Console.Write(accept_Q_times_request_ALL_total_control[all_control] + " =  accept_Q_times_request_ALL_total_control  ");
            Console.WriteLine();
            Console.Write(accept_Q_times_request_FT_total_control[all_control] + " =   accept_Q_times_request_FT_total_control ");
            Console.WriteLine();
            Console.Write(Q_times_request_ALL_total_control[all_control] + "  =  Q_times_request_ALL_total_control ");
            Console.WriteLine();
            Console.Write(RATE_accept_Q_times_request_pickup_FT_total_control[all_control] + " =   RATE_accept_Q_times_request_pickup_FT_total_control ");
            Console.WriteLine();
            Console.Write(RATE_accept_Q_times_request_dropoff_FT_total_control[all_control] + "  =   RATE_accept_Q_times_request_dropoff_FT_total_control");
            Console.WriteLine();
            Console.Write(RATE_accept_Q_times_request_ALL_total_control[all_control] + " =  RATE_accept_Q_times_request_ALL_total_control  ");
            Console.WriteLine();
            Console.Write(RATE_accept_Q_times_request_FT_total_control[all_control] + "  =  RATE_accept_Q_times_request_FT_total_control ");
            Console.WriteLine();
        }
 
        public void Floating_target_results_Excel_Table1029_ALL0_220419(int TYPE)
        {

            input_data_Manhattan_Intersection_Coordinates_Mod_Dis_Car();

            Console.WriteLine("OK!");
            Console.WriteLine("ALL 14 days");

            double AVG_FT_vehicle_miles_driving_control = 0;
            double AVG_FT_vehicle_miles_peracceptedcustomer_driving_control = 0;

            Console.WriteLine("TYPE = " + TYPE);
            for (int all_control = 1; all_control < 14 + 1; all_control++)
            {
                Floating_target_results_Excel_Table1029_220419(all_control, TYPE);
                Floating_target_results_Excel_Table1_220322_220419(all_control, TYPE);
            }

            for (int all_control = 1; all_control < 14 + 1; all_control++)
            {
                Console.Write(home_vehicle_miles_driving_control[all_control] * 0.000621371192 + "    ");

            }

            Console.WriteLine();

            for (int all_control = 1; all_control < 14 + 1; all_control++)
            {
                Console.Write(FT_vehicle_miles_driving_control[all_control] * 0.000621371192 + "    ");

                AVG_FT_vehicle_miles_driving_control = AVG_FT_vehicle_miles_driving_control + FT_vehicle_miles_driving_control[all_control] * 0.000621371192;
                AVG_FT_vehicle_miles_peracceptedcustomer_driving_control = AVG_FT_vehicle_miles_peracceptedcustomer_driving_control + FT_vehicle_miles_driving_control[all_control] * 0.000621371192 / Acceptance_control[all_control];
            }
            Console.WriteLine();

            AVG_FT_vehicle_miles_driving_control = AVG_FT_vehicle_miles_driving_control / 14;
            AVG_FT_vehicle_miles_peracceptedcustomer_driving_control = AVG_FT_vehicle_miles_peracceptedcustomer_driving_control / 14;

            Console.WriteLine("AVG_FT_vehicle_miles_driving_control = " + AVG_FT_vehicle_miles_driving_control);
            Console.WriteLine("AVG_FT_vehicle_miles_peracceptedcustomer_driving_control = " + AVG_FT_vehicle_miles_peracceptedcustomer_driving_control);

        }

        public void input_all_fare_0921(int date)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            if (date == 1)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_01_2019_Mod.mdb";
            if (date == 2)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_02_2019_Mod.mdb";
            if (date == 3)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_03_2019_Mod.mdb";
            if (date == 4)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_04_2019_Mod.mdb";
            if (date == 5)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_05_2019_Mod.mdb";
            if (date == 6)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_06_2019_Mod.mdb";
            if (date == 7)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_07_2019_Mod.mdb";
            if (date == 8)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_08_2019_Mod.mdb";
            if (date == 9)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_09_2019_Mod.mdb";
            if (date == 10)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_10_2019_Mod.mdb";
            if (date == 11)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_11_2019_Mod.mdb";
            if (date == 12)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_12_2019_Mod.mdb";
            if (date == 13)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_13_2019_Mod.mdb";
            if (date == 14)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_14_2019_Mod.mdb";
            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Instance_file";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();

            int jj = 0;
            while (myReader001.Read())
            {
                if (myReader001.GetDouble(2) <= 1.5 * 3600 && myReader001.GetDouble(2) > 0)
                {
                    jj = myReader001.GetInt16(1);
                    g_truej[jj] = myReader001.GetDouble(4);
                }

                if (myReader001.GetDouble(2) > 1.5 * 3600)
                    break;
            }
            myReader001.Close();
            objConnection.Close();
        }

        public void input_all_Q_0921_220422(int date)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
            if (date == 1)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_01_2019_Mod.mdb";
            if (date == 2)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_02_2019_Mod.mdb";
            if (date == 3)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_03_2019_Mod.mdb";
            if (date == 4)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_04_2019_Mod.mdb";
            if (date == 5)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_05_2019_Mod.mdb";
            if (date == 6)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_06_2019_Mod.mdb";
            if (date == 7)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_07_2019_Mod.mdb";
            if (date == 8)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_08_2019_Mod.mdb";
            if (date == 9)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_09_2019_Mod.mdb";
            if (date == 10)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_10_2019_Mod.mdb";
            if (date == 11)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_11_2019_Mod.mdb";
            if (date == 12)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_12_2019_Mod.mdb";
            if (date == 13)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_13_2019_Mod.mdb";
            if (date == 14)
                strConnection += @"Data Source=C:\Users\wzhan\Desktop\MScomments\online_routing\Instances_NYC_taxi_Mod\Instances_NYC_taxi_Mod\Instance_file_12_14_2019_Mod.mdb";
            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            string str1 = "select * from Instance_file";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();

            int jj = 0;
            while (myReader001.Read())
            {
                if (myReader001.GetDouble(2) <= 1.5 * 3600 && myReader001.GetDouble(2) > 0)
                {
                    jj = myReader001.GetInt16(1);
                    q_truej[jj] = myReader001.GetInt16(3);
                }

                if (myReader001.GetDouble(2) > 1.5 * 3600)
                    break;
            }
            myReader001.Close();
            objConnection.Close();
        }
       
        public void Floating_target_results_Excel_Table1_ALL_220322_220419(int TYPE)
        {
            Console.WriteLine("TYPE = " + TYPE);
            double AVG_Acceptance_control = 0;
            double AVG_Acceptance_rate_control = 0;
            double AVG_FT_Pickup_rate_control = 0;
            double AVG_FT_Dropoff_rate_control = 0;
            double AVG_Shared_path_rate_control = 0;


            for (int all_control = 1; all_control < 14 + 1; all_control++)
            {
                {
                    Floating_target_results_Excel_Table1_220322_220419(all_control, TYPE); 
                    Console.Write(Total_customer_control[all_control] + "    " + Acceptance_control[all_control] + "    " + Acceptance_rate_control[all_control] + "    " +
                       All_Floating_control[all_control] + "    " + FT_Pickup_control[all_control] + "    " + FT_Dropoff_control[all_control] + "    " +
                       All_Floating_rate_control[all_control] + "    " + FT_Pickup_rate_control[all_control] + "    " + FT_Dropoff_rate_control[all_control] + "    " +
                       Total_unreposition_path_control[all_control] + "    " + Shared_path_control[all_control] + "    " + Shared_path_rate_control[all_control]);
                    Console.WriteLine();
                }
                AVG_Acceptance_control = AVG_Acceptance_control + Acceptance_control[all_control];
                AVG_Acceptance_rate_control = AVG_Acceptance_rate_control + Acceptance_rate_control[all_control];
                AVG_FT_Pickup_rate_control = AVG_FT_Pickup_rate_control + FT_Pickup_rate_control[all_control];
                AVG_FT_Dropoff_rate_control = AVG_FT_Dropoff_rate_control + FT_Dropoff_rate_control[all_control];
                AVG_Shared_path_rate_control = AVG_Shared_path_rate_control + Shared_path_rate_control[all_control];

             
            }
            AVG_Acceptance_control = AVG_Acceptance_control / 14;
            AVG_Acceptance_rate_control = AVG_Acceptance_rate_control / 14;
            AVG_FT_Pickup_rate_control = AVG_FT_Pickup_rate_control / 14;
            AVG_FT_Dropoff_rate_control = AVG_FT_Dropoff_rate_control / 14;
            AVG_Shared_path_rate_control = AVG_Shared_path_rate_control / 14;

        
            Console.WriteLine("AVG_Acceptance_control = " + AVG_Acceptance_control);
            Console.WriteLine("AVG_Acceptance_rate_control = " + AVG_Acceptance_rate_control);
            Console.WriteLine("AVG_FT_Pickup_rate_control = " + AVG_FT_Pickup_rate_control);
            Console.WriteLine("AVG_FT_Dropoff_rate_control = " + AVG_FT_Dropoff_rate_control);
            Console.WriteLine("AVG_Shared_path_rate_control = " + AVG_Shared_path_rate_control);

        }

        public void Floating_target_results_Excel_Table0524_ALL_220419(int TYPE)
        {
            Console.WriteLine("TYPE = " + TYPE);
            double AVG_Total_profit_control = 0;
            double AVG_Total_fare_control = 0;
            double AVG_Total_cost_control = 0;
            double AVG_average_cost_peracceptedcustomer_control = 0;

            for (int all_control = 1; all_control < 14 + 1; all_control++)
            {
                Floating_target_results_Excel_Table_0524_220419(all_control, TYPE);

                Console.Write(Total_profit_control[all_control] + "    " + Total_fare_control[all_control] + "    " + Total_cost_customer_control[all_control] + "    " + Total_cost_reposition_control[all_control]);
                Console.WriteLine();
                AVG_Total_profit_control = AVG_Total_profit_control + Total_profit_control[all_control];
                AVG_Total_fare_control = AVG_Total_fare_control + Total_fare_control[all_control];

                Floating_target_results_Excel_Table1_220322_220419(all_control, TYPE);
                Console.WriteLine((Total_fare_control[all_control] - Total_profit_control[all_control]) / Acceptance_control[all_control]);
                AVG_average_cost_peracceptedcustomer_control = AVG_average_cost_peracceptedcustomer_control + (Total_fare_control[all_control] - Total_profit_control[all_control]) / Acceptance_control[all_control];

            }

            AVG_Total_profit_control = AVG_Total_profit_control / 14;
            AVG_Total_fare_control = AVG_Total_fare_control / 14;
            AVG_Total_cost_control = AVG_Total_fare_control - AVG_Total_profit_control;
            Console.WriteLine("AVG_Total_profit_control = " + AVG_Total_profit_control);
            Console.WriteLine("AVG_Total_fare_control = " + AVG_Total_fare_control);
            Console.WriteLine("AVG_Total_cost_control = " + AVG_Total_cost_control);

            AVG_average_cost_peracceptedcustomer_control = AVG_average_cost_peracceptedcustomer_control / 14;
            Console.WriteLine("AVG_average_cost_peracceptedcustomer_control = " + AVG_average_cost_peracceptedcustomer_control);

        }

  
        public void Floating_target_results_Excel_Table0524_ALL_220419_220420_220421()
        {
            double AVG_Profit_increase9 = 0;
            double AVG_Revenue_increase9 = 0;
            double AVG_Cost_decrease9 = 0;

            double AVG_Profit_increase5 = 0;
            double AVG_Revenue_increase5 = 0;
            double AVG_Cost_decrease5 = 0;

            double AVG_Profit_increase6 = 0;
            double AVG_Revenue_increase6 = 0;
            double AVG_Cost_decrease6 = 0;

            double AVG_Profit_increase4 = 0;
            double AVG_Revenue_increase4 = 0;
            double AVG_Cost_decrease4 = 0;

            double AVG_downstream_cost_savings9 = 0;
            double AVG_upstream_change_in_customer_mix9 = 0;

            double AVG_downstream_cost_savings5 = 0;
            double AVG_upstream_change_in_customer_mix5 = 0;

            double AVG_downstream_cost_savings6 = 0;
            double AVG_upstream_change_in_customer_mix6 = 0;

            double AVG_downstream_cost_savings4 = 0;
            double AVG_upstream_change_in_customer_mix4 = 0;

            double AVG_max_discount_uniform9 = 0;
            double AVG_max_discount_uniform5 = 0;
            double AVG_max_discount_uniform6 = 0;
            double AVG_max_discount_uniform4 = 0;

            double AVG_max_discount_targeted9 = 0;
            double AVG_max_discount_targeted5 = 0;
            double AVG_max_discount_targeted6 = 0;
            double AVG_max_discount_targeted4 = 0;


            double AVG_max_discount_prorated9 = 0;
            double AVG_max_discount_prorated5 = 0;
            double AVG_max_discount_prorated6 = 0;
            double AVG_max_discount_prorated4 = 0;


            double AVG_RATE_extra_customers_withQ_4_v1 = 0;
            double AVG_RATE_extra_customers_withQ_4_v2 = 0;

            double AVG_RATE_accept_Q_times_request_pickup_FT_total = 0;
            double AVG_RATE_accept_Q_times_request_dropoff_FT_total = 0;

            double AVG_walking_distance4 = 0;

            for (int TYPE = 1201; TYPE < 1214 + 1; TYPE++)
            {
                Console.WriteLine("TYPE = " + TYPE);
                for (int all_control = 1; all_control < 16 + 1; all_control++)
                {

                    if (all_control == 9)
                    {
                        Floating_target_results_Excel_Table_0524_220419(9, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(9, TYPE);
    
                        Floating_target_results_Excel_Table_0524_220419(1, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(1, TYPE);
            

                        AVG_Profit_increase9 = AVG_Profit_increase9 + (Total_profit_control[9] - Total_profit_control[1]) / Total_profit_control[1];
                        AVG_Revenue_increase9 = AVG_Revenue_increase9 + (Total_fare_control[9] - Total_fare_control[1]) / Total_fare_control[1];
                        AVG_Cost_decrease9 = AVG_Cost_decrease9 + (Total_cost_customer_control[9] + Total_cost_reposition_control[9] - Total_cost_customer_control[1] - Total_cost_reposition_control[1]) / (Total_cost_customer_control[1] + Total_cost_reposition_control[1]);

                        double C0 = (Total_cost_customer_control[1] + Total_cost_reposition_control[1]) / Acceptance_control[1];
                        double C1 = (Total_cost_customer_control[9] + Total_cost_reposition_control[9]) / Acceptance_control[9];
                        double Q0 = Acceptance_control[1];
                        double Q1 = Acceptance_control[9];
                        double P0 = Total_fare_control[1] / Acceptance_control[1];
                        double P1 = Total_fare_control[9] / Acceptance_control[9];
                        double Pai_D = (C0 - C1) * Q0;
                        double Pai_U = (P1 * Q1 - P0 * Q0) - C1 * (Q1 - Q0);
                        double Pai = Pai_D + Pai_U;

                        Console.WriteLine("CHECK: " + Pai + " = " + (Total_profit_control[9] - Total_profit_control[1]));
                        AVG_downstream_cost_savings9 = AVG_downstream_cost_savings9 + Pai_D / Pai;
                        AVG_upstream_change_in_customer_mix9 = AVG_upstream_change_in_customer_mix9 + Pai_U / Pai;
                        Console.WriteLine(" Pai_D / Pai = " + Pai_D / Pai);
                        Console.WriteLine(" Pai_U / Pai = " + Pai_U / Pai);

                        //================================================
                        double Discount_uniform_thisone = 1 - (Total_profit_control[1] + (Total_cost_customer_control[9] + Total_cost_reposition_control[9])) / Total_fare_control[9];
                        AVG_max_discount_uniform9 = AVG_max_discount_uniform9 + Discount_uniform_thisone;

                        //================================================
                        Floating_target_results_Excel_Table0921_ALL_220421(9, TYPE);
                        Floating_target_results_Excel_Table0921_ALL_220421(1, TYPE);

                        double Discount_targeted_thisone = 1 - (Total_profit_control[1] + (Total_cost_customer_control[9] + Total_cost_reposition_control[9]) - (Total_fare_control[9] - Fare_FT_control[9])) / Fare_FT_control[9];
                        AVG_max_discount_targeted9 = AVG_max_discount_targeted9 + Discount_targeted_thisone;

                        double Discount_prorated_thisone = (Total_profit_control[9] - Total_profit_control[1]) / (pickup_WD_total_control[9] + dropoff_WD_total_control[9]);
                        AVG_max_discount_prorated9 = AVG_max_discount_prorated9 + Discount_prorated_thisone;
                    }
                    if (all_control == 5)
                    {
                        Floating_target_results_Excel_Table_0524_220419(5, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(5, TYPE);
       
                        Floating_target_results_Excel_Table_0524_220419(3, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(3, TYPE);
           
                        AVG_Profit_increase5 = AVG_Profit_increase5 + (Total_profit_control[5] - Total_profit_control[3]) / Total_profit_control[3];
                        AVG_Revenue_increase5 = AVG_Revenue_increase5 + (Total_fare_control[5] - Total_fare_control[3]) / Total_fare_control[3];
                        AVG_Cost_decrease5 = AVG_Cost_decrease5 + (Total_cost_customer_control[5] + Total_cost_reposition_control[5] - Total_cost_customer_control[3] - Total_cost_reposition_control[3]) / (Total_cost_customer_control[3] + Total_cost_reposition_control[3]);

                        double C0 = (Total_cost_customer_control[3] + Total_cost_reposition_control[3]) / Acceptance_control[3];
                        double C1 = (Total_cost_customer_control[5] + Total_cost_reposition_control[5]) / Acceptance_control[5];
                        double Q0 = Acceptance_control[3];
                        double Q1 = Acceptance_control[5];
                        double P0 = Total_fare_control[3] / Acceptance_control[3];
                        double P1 = Total_fare_control[5] / Acceptance_control[5];
                        double Pai_D = (C0 - C1) * Q0;
                        double Pai_U = (P1 * Q1 - P0 * Q0) - C1 * (Q1 - Q0);
                        double Pai = Pai_D + Pai_U;

                        Console.WriteLine("CHECK: " + Pai + " = " + (Total_profit_control[5] - Total_profit_control[3]));
                        AVG_downstream_cost_savings5 = AVG_downstream_cost_savings5 + Pai_D / Pai;
                        AVG_upstream_change_in_customer_mix5 = AVG_upstream_change_in_customer_mix5 + Pai_U / Pai;
                        Console.WriteLine(" Pai_D / Pai = " + Pai_D / Pai);
                        Console.WriteLine(" Pai_U / Pai = " + Pai_U / Pai);

                        double Discount_uniform_thisone = 1 - (Total_profit_control[3] + (Total_cost_customer_control[5] + Total_cost_reposition_control[5])) / Total_fare_control[5];
                        AVG_max_discount_uniform5 = AVG_max_discount_uniform5 + Discount_uniform_thisone;

                        //================================================
                        Floating_target_results_Excel_Table0921_ALL_220421(5, TYPE);
                        Floating_target_results_Excel_Table0921_ALL_220421(3, TYPE);

                        double Discount_targeted_thisone = 1 - (Total_profit_control[3] + (Total_cost_customer_control[5] + Total_cost_reposition_control[5]) - (Total_fare_control[5] - Fare_FT_control[5])) / Fare_FT_control[5];
                        AVG_max_discount_targeted5 = AVG_max_discount_targeted5 + Discount_targeted_thisone;

                        double Discount_prorated_thisone = (Total_profit_control[5] - Total_profit_control[3]) / (pickup_WD_total_control[5] + dropoff_WD_total_control[5]);
                        AVG_max_discount_prorated5 = AVG_max_discount_prorated5 + Discount_prorated_thisone;

                    }
                    if (all_control == 6)
                    {
                        Floating_target_results_Excel_Table_0524_220419(6, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(6, TYPE);
            
                        Floating_target_results_Excel_Table_0524_220419(3, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(3, TYPE);
       

                        AVG_Profit_increase6 = AVG_Profit_increase6 + (Total_profit_control[6] - Total_profit_control[3]) / Total_profit_control[3];
                        AVG_Revenue_increase6 = AVG_Revenue_increase6 + (Total_fare_control[6] - Total_fare_control[3]) / Total_fare_control[3];
                        AVG_Cost_decrease6 = AVG_Cost_decrease6 + (Total_cost_customer_control[6] + Total_cost_reposition_control[6] - Total_cost_customer_control[3] - Total_cost_reposition_control[3]) / (Total_cost_customer_control[3] + Total_cost_reposition_control[3]);



                        double C0 = (Total_cost_customer_control[3] + Total_cost_reposition_control[3]) / Acceptance_control[3];
                        double C1 = (Total_cost_customer_control[6] + Total_cost_reposition_control[6]) / Acceptance_control[6];
                        double Q0 = Acceptance_control[3];
                        double Q1 = Acceptance_control[6];
                        double P0 = Total_fare_control[3] / Acceptance_control[3];
                        double P1 = Total_fare_control[6] / Acceptance_control[6];
                        double Pai_D = (C0 - C1) * Q0;
                        double Pai_U = (P1 * Q1 - P0 * Q0) - C1 * (Q1 - Q0);
                        double Pai = Pai_D + Pai_U;

                        Console.WriteLine("CHECK: " + Pai + " = " + (Total_profit_control[6] - Total_profit_control[3]));
                        AVG_downstream_cost_savings6 = AVG_downstream_cost_savings6 + Pai_D / Pai;
                        AVG_upstream_change_in_customer_mix6 = AVG_upstream_change_in_customer_mix6 + Pai_U / Pai;
                        Console.WriteLine(" Pai_D / Pai = " + Pai_D / Pai);
                        Console.WriteLine(" Pai_U / Pai = " + Pai_U / Pai);

                        double Discount_uniform_thisone = 1 - (Total_profit_control[3] + (Total_cost_customer_control[6] + Total_cost_reposition_control[6])) / Total_fare_control[6];
                        AVG_max_discount_uniform6 = AVG_max_discount_uniform6 + Discount_uniform_thisone;

                        //================================================
                        Floating_target_results_Excel_Table0921_ALL_220421(6, TYPE);
                        Floating_target_results_Excel_Table0921_ALL_220421(3, TYPE);

                        double Discount_targeted_thisone = 1 - (Total_profit_control[3] + (Total_cost_customer_control[6] + Total_cost_reposition_control[6]) - (Total_fare_control[6] - Fare_FT_control[6])) / Fare_FT_control[6];
                        AVG_max_discount_targeted6 = AVG_max_discount_targeted6 + Discount_targeted_thisone;

                        double Discount_prorated_thisone = (Total_profit_control[6] - Total_profit_control[3]) / (pickup_WD_total_control[6] + dropoff_WD_total_control[6]);
                        AVG_max_discount_prorated6 = AVG_max_discount_prorated6 + Discount_prorated_thisone;

                    }

                    if (all_control == 4)
                    {
                        Floating_target_results_Excel_Table_0524_220419(4, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(4, TYPE);

                        Floating_target_results_Excel_Table_0524_220419(3, TYPE);
                        Floating_target_results_Excel_Table1_220322_220419(3, TYPE);
               
                        AVG_Profit_increase4 = AVG_Profit_increase4 + (Total_profit_control[4] - Total_profit_control[3]) / Total_profit_control[3];
                        AVG_Revenue_increase4 = AVG_Revenue_increase4 + (Total_fare_control[4] - Total_fare_control[3]) / Total_fare_control[3];
                        AVG_Cost_decrease4 = AVG_Cost_decrease4 + (Total_cost_customer_control[4] + Total_cost_reposition_control[4] - Total_cost_customer_control[3] - Total_cost_reposition_control[3]) / (Total_cost_customer_control[3] + Total_cost_reposition_control[3]);


                        double C0 = (Total_cost_customer_control[3] + Total_cost_reposition_control[3]) / Acceptance_control[3];
                        double C1 = (Total_cost_customer_control[4] + Total_cost_reposition_control[4]) / Acceptance_control[4];
                        double Q0 = Acceptance_control[3];
                        double Q1 = Acceptance_control[4];
                        double P0 = Total_fare_control[3] / Acceptance_control[3];
                        double P1 = Total_fare_control[4] / Acceptance_control[4];
                        double Pai_D = (C0 - C1) * Q0;
                        double Pai_U = (P1 * Q1 - P0 * Q0) - C1 * (Q1 - Q0);
                        double Pai = Pai_D + Pai_U;

                        Console.WriteLine("CHECK: " + Pai + " = " + (Total_profit_control[4] - Total_profit_control[3]));
                        AVG_downstream_cost_savings4 = AVG_downstream_cost_savings4 + Pai_D / Pai;
                        AVG_upstream_change_in_customer_mix4 = AVG_upstream_change_in_customer_mix4 + Pai_U / Pai;
                        Console.WriteLine(" Pai_D / Pai = " + Pai_D / Pai);
                        Console.WriteLine(" Pai_U / Pai = " + Pai_U / Pai);


                        double Discount_uniform_thisone = 1 - (Total_profit_control[3] + (Total_cost_customer_control[4] + Total_cost_reposition_control[4])) / Total_fare_control[4];
                        AVG_max_discount_uniform4 = AVG_max_discount_uniform4 + Discount_uniform_thisone;

                        //================================================
                        Floating_target_results_Excel_Table0921_ALL_220421(4, TYPE);
                        Floating_target_results_Excel_Table0921_ALL_220421(3, TYPE);

                        double Discount_targeted_thisone = 1 - (Total_profit_control[3] + (Total_cost_customer_control[4] + Total_cost_reposition_control[4]) - (Total_fare_control[4] - Fare_FT_control[4])) / Fare_FT_control[4];
                        AVG_max_discount_targeted4 = AVG_max_discount_targeted4 + Discount_targeted_thisone;

                        double Discount_prorated_thisone = (Total_profit_control[4] - Total_profit_control[3]) / (pickup_WD_total_control[4] + dropoff_WD_total_control[4]);
                        AVG_max_discount_prorated4 = AVG_max_discount_prorated4 + Discount_prorated_thisone;

                        //=================================================================================
                        Floating_target_results_Excel_Table1_220322_220419(4, TYPE);
                        AVG_walking_distance4 = AVG_walking_distance4 + pickup_WD_total_control[4] / FT_Pickup_control[4] + dropoff_WD_total_control[4] / FT_Dropoff_control[4];


                        AVG_RATE_extra_customers_withQ_4_v1 = AVG_RATE_extra_customers_withQ_4_v1 + (RATE_accept_Q_times_request_ALL_total_control[4] - RATE_accept_Q_times_request_ALL_total_control[3]);
                        AVG_RATE_extra_customers_withQ_4_v2 = AVG_RATE_extra_customers_withQ_4_v2 + (accept_Q_times_request_ALL_total_control[4] - accept_Q_times_request_ALL_total_control[3]) / accept_Q_times_request_ALL_total_control[3];
                        AVG_RATE_accept_Q_times_request_pickup_FT_total = AVG_RATE_accept_Q_times_request_pickup_FT_total + RATE_accept_Q_times_request_pickup_FT_total_control[4];
                        AVG_RATE_accept_Q_times_request_dropoff_FT_total = AVG_RATE_accept_Q_times_request_dropoff_FT_total + RATE_accept_Q_times_request_dropoff_FT_total_control[4];
                    }
                }
            }
            AVG_Profit_increase9 = AVG_Profit_increase9 / 14;
            AVG_Revenue_increase9 = AVG_Revenue_increase9 / 14;
            AVG_Cost_decrease9 = AVG_Cost_decrease9 / 14;
            Console.WriteLine("AVG_Profit_increase9 = " + AVG_Profit_increase9);
            Console.WriteLine("AVG_Revenue_increase9 = " + AVG_Revenue_increase9);
            Console.WriteLine("AVG_Cost_decrease9 = " + -AVG_Cost_decrease9);

            AVG_Profit_increase5 = AVG_Profit_increase5 / 14;
            AVG_Revenue_increase5 = AVG_Revenue_increase5 / 14;
            AVG_Cost_decrease5 = AVG_Cost_decrease5 / 14;
            Console.WriteLine("AVG_Profit_increase5 = " + AVG_Profit_increase5);
            Console.WriteLine("AVG_Revenue_increase5 = " + AVG_Revenue_increase5);
            Console.WriteLine("AVG_Cost_decrease5 = " + -AVG_Cost_decrease5);


            AVG_Profit_increase6 = AVG_Profit_increase6 / 14;
            AVG_Revenue_increase6 = AVG_Revenue_increase6 / 14;
            AVG_Cost_decrease6 = AVG_Cost_decrease6 / 14;
            Console.WriteLine("AVG_Profit_increase6 = " + AVG_Profit_increase6);
            Console.WriteLine("AVG_Revenue_increase6 = " + AVG_Revenue_increase6);
            Console.WriteLine("AVG_Cost_decrease6 = " + -AVG_Cost_decrease6);


            AVG_Profit_increase4 = AVG_Profit_increase4 / 14;
            AVG_Revenue_increase4 = AVG_Revenue_increase4 / 14;
            AVG_Cost_decrease4 = AVG_Cost_decrease4 / 14;
            Console.WriteLine("AVG_Profit_increase4 = " + AVG_Profit_increase4);
            Console.WriteLine("AVG_Revenue_increase4 = " + AVG_Revenue_increase4);
            Console.WriteLine("AVG_Cost_decrease4 = " + -AVG_Cost_decrease4);



            AVG_downstream_cost_savings9 = AVG_downstream_cost_savings9 / 14;
            AVG_upstream_change_in_customer_mix9 = AVG_upstream_change_in_customer_mix9 / 14;
            Console.WriteLine("AVG_downstream_cost_savings9 = " + AVG_downstream_cost_savings9);
            Console.WriteLine("AVG_upstream_change_in_customer_mix9 = " + AVG_upstream_change_in_customer_mix9);

            AVG_downstream_cost_savings5 = AVG_downstream_cost_savings5 / 14;
            AVG_upstream_change_in_customer_mix5 = AVG_upstream_change_in_customer_mix5 / 14;
            Console.WriteLine("AVG_downstream_cost_savings5 = " + AVG_downstream_cost_savings5);
            Console.WriteLine("AVG_upstream_change_in_customer_mix5 = " + AVG_upstream_change_in_customer_mix5);

            AVG_downstream_cost_savings6 = AVG_downstream_cost_savings6 / 14;
            AVG_upstream_change_in_customer_mix6 = AVG_upstream_change_in_customer_mix6 / 14;
            Console.WriteLine("AVG_downstream_cost_savings6 = " + AVG_downstream_cost_savings6);
            Console.WriteLine("AVG_upstream_change_in_customer_mix6 = " + AVG_upstream_change_in_customer_mix6);

            AVG_downstream_cost_savings4 = AVG_downstream_cost_savings4 / 14;
            AVG_upstream_change_in_customer_mix4 = AVG_upstream_change_in_customer_mix4 / 14;
            Console.WriteLine("AVG_downstream_cost_savings4 = " + AVG_downstream_cost_savings4);
            Console.WriteLine("AVG_upstream_change_in_customer_mix4 = " + AVG_upstream_change_in_customer_mix4);

            //===========================================================================================
            AVG_max_discount_uniform9 = AVG_max_discount_uniform9 / 14;
            AVG_max_discount_targeted9 = AVG_max_discount_targeted9 / 14;
            AVG_max_discount_prorated9 = AVG_max_discount_prorated9 / 14;
            Console.WriteLine("AVG_max_discount_uniform9 = " + AVG_max_discount_uniform9);
            Console.WriteLine("AVG_max_discount_targeted9 = " + AVG_max_discount_targeted9);
            Console.WriteLine("AVG_max_discount_prorated9 = " + AVG_max_discount_prorated9);

            AVG_max_discount_uniform5 = AVG_max_discount_uniform5 / 14;
            AVG_max_discount_targeted5 = AVG_max_discount_targeted5 / 14;
            AVG_max_discount_prorated5 = AVG_max_discount_prorated5 / 14;
            Console.WriteLine("AVG_max_discount_uniform5 = " + AVG_max_discount_uniform5);
            Console.WriteLine("AVG_max_discount_targeted5 = " + AVG_max_discount_targeted5);
            Console.WriteLine("AVG_max_discount_prorated5 = " + AVG_max_discount_prorated5);

            AVG_max_discount_uniform6 = AVG_max_discount_uniform6 / 14;
            AVG_max_discount_targeted6 = AVG_max_discount_targeted6 / 14;
            AVG_max_discount_prorated6 = AVG_max_discount_prorated6 / 14;
            Console.WriteLine("AVG_max_discount_uniform6 = " + AVG_max_discount_uniform6);
            Console.WriteLine("AVG_max_discount_targeted6 = " + AVG_max_discount_targeted6);
            Console.WriteLine("AVG_max_discount_prorated6 = " + AVG_max_discount_prorated6);

            AVG_max_discount_uniform4 = AVG_max_discount_uniform4 / 14;
            AVG_max_discount_targeted4 = AVG_max_discount_targeted4 / 14;
            AVG_max_discount_prorated4 = AVG_max_discount_prorated4 / 14;
            Console.WriteLine("AVG_max_discount_uniform4 = " + AVG_max_discount_uniform4);
            Console.WriteLine("AVG_max_discount_targeted4 = " + AVG_max_discount_targeted4);
            Console.WriteLine("AVG_max_discount_prorated4 = " + AVG_max_discount_prorated4);

        }


        public void Floating_target_results_Excel_Table_0524_220419(int all_control, int TYPE)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";

     
            {
                if (TYPE == 1)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                }
                if (TYPE == 9)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                }
                if (TYPE == 3)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                }
                if (TYPE == 5)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                }
                if (TYPE == 6)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                }
                if (TYPE == 4)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }

                if (TYPE == 1201)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1202)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1203)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1204)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1205)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1206)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1207)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1208)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1209)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1210)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1211)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1212)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1213)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1214)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }
            }

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            //================================================================================================


            int Total_subpath = 0;
            int Total_reposition_subpath = 0;
            int Total_noreposition_subpath = 0;
            int Total_share_subpath = 0;
            int Total_unshare_subpath = 0;
            int[] Share_subpath = new int[100000];

            double[] Fare_subpath = new double[100000];
            double[] Cost_customer_subpath = new double[100000];
            double[] Cost_reposition_subpath = new double[100000];

            double Total_Duration_not_shared_subpath = 0;
            double Total_Duration_shared_subpath = 0;
            double Total_Duration_all_subpath = 0;

            double Total_fare = 0;
            double Total_cost_customer = 0;
            double Total_cost_reposition = 0;

            string str3 = "select * from subpath_profits";
            OleDbCommand myCommand003 = new OleDbCommand(str3, objConnection);
            OleDbDataReader myReader003 = myCommand003.ExecuteReader();

            while (myReader003.Read())
            {
                int dd3 = myReader003.GetInt32(1);
                Share_subpath[dd3] = myReader003.GetInt16(10); 

                Fare_subpath[dd3] = myReader003.GetDouble(7);
                if (Fare_subpath[dd3] > 0)
                {
                    Cost_customer_subpath[dd3] = myReader003.GetDouble(8);
                }
                else
                {
                    Cost_reposition_subpath[dd3] = myReader003.GetDouble(8);
                }

                if (dd3 > Total_subpath)
                    Total_subpath = dd3;
            }
            myReader003.Close();

            for (int ppp = 1; ppp < Total_subpath + 1; ppp++)
            {
                if (Share_subpath[ppp] == 0)
                {
                    Total_reposition_subpath = Total_reposition_subpath + 1;
                }
                else
                {
                    Total_noreposition_subpath = Total_noreposition_subpath + 1;
                    Total_Duration_all_subpath = Total_Duration_all_subpath + ((double)Cost_customer_subpath[ppp] / 3.2);

                    if (Share_subpath[ppp] == 1)
                    {
                        Total_unshare_subpath = Total_unshare_subpath + 1;
                        Total_Duration_not_shared_subpath = Total_Duration_not_shared_subpath + ((double)Cost_customer_subpath[ppp] / 3.2);
                    }
                    if (Share_subpath[ppp] > 1)
                    {
                        Total_share_subpath = Total_share_subpath + 1;
                        Total_Duration_shared_subpath = Total_Duration_shared_subpath + ((double)Cost_customer_subpath[ppp] / 3.2);
                    }
                }

                Total_fare = Total_fare + Fare_subpath[ppp];
                Total_cost_customer = Total_cost_customer + Cost_customer_subpath[ppp];
                Total_cost_reposition = Total_cost_reposition + Cost_reposition_subpath[ppp];
            }

            Total_unreposition_path_control[all_control] = Total_noreposition_subpath;
            Shared_path_control[all_control] = Total_share_subpath;
            Shared_path_rate_control[all_control] = ((double)Total_share_subpath / (double)Total_noreposition_subpath);


            Total_fare_control[all_control] = Total_fare;
            Total_cost_customer_control[all_control] = Total_cost_customer;
            Total_cost_reposition_control[all_control] = Total_cost_reposition;
            Total_profit_control[all_control] = Total_fare - Total_cost_customer - Total_cost_reposition;

            //================================================================================================
            Duration_all_subpath_control[all_control] = (double)Total_Duration_all_subpath / (double)Total_noreposition_subpath;
            Duration_not_shared_subpath_control[all_control] = (double)Total_Duration_not_shared_subpath / ((double)Total_noreposition_subpath - (double)Total_share_subpath);
            Duration_shared_subpath_control[all_control] = (double)Total_Duration_shared_subpath / (double)Total_share_subpath;

            objConnection.Close();

        }

        public double[] pickup_WD_total_control = new double[20];
        public double[] dropoff_WD_total_control = new double[20];
        public double[] Fare_FT_control = new double[20];

        public double[] accept_Q_times_request_pickup_FT_total_control = new double[20];
        public double[] accept_Q_times_request_dropoff_FT_total_control = new double[20];
        public double[] accept_Q_times_request_ALL_total_control = new double[20];
        public double[] accept_Q_times_request_FT_total_control = new double[20];
        public double[] Q_times_request_ALL_total_control = new double[20];

        public double[] RATE_accept_Q_times_request_pickup_FT_total_control = new double[20];
        public double[] RATE_accept_Q_times_request_dropoff_FT_total_control = new double[20];
        public double[] RATE_accept_Q_times_request_FT_total_control = new double[20];
        public double[] RATE_accept_Q_times_request_ALL_total_control = new double[20];

        public double[] home_vehicle_miles_driving_control = new double[20];
        public double[] FT_vehicle_miles_driving_control = new double[20];


        public double[] average_delay_control = new double[20];
        public int[] Averge_delay_number_control = new int[20];

        public void Floating_target_results_Excel_Table0921_220421(int all_control, int TYPE)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
        
            {
                if (TYPE == 1)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                }
                if (TYPE == 9)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                }
                if (TYPE == 3)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                }
                if (TYPE == 5)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                }
                if (TYPE == 6)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                }
                if (TYPE == 4)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }

                if (TYPE == 1201)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1202)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1203)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1204)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1205)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1206)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1207)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1208)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1209)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1210)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1211)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1212)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1213)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1214)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }
            }

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();


            int[] Status_customer = new int[20000];
            int Total_customer = 0;
            int Total_accept_customer = 0;

            double[] pickup_WD = new double[20000];
            double[] dropoff_WD = new double[20000];
            double pickup_WD_total = 0;
            double dropoff_WD_total = 0;
            double Fare_FT = 0;


            string str2 = "select * from customer_status";
            OleDbCommand myCommand002 = new OleDbCommand(str2, objConnection);
            OleDbDataReader myReader002 = myCommand002.ExecuteReader();

            while (myReader002.Read())
            {
                int dd2 = myReader002.GetInt32(1);
                Status_customer[dd2] = myReader002.GetInt16(2);

                pickup_WD[dd2] = myReader002.GetDouble(6);
                dropoff_WD[dd2] = myReader002.GetDouble(14);

                if (dd2 > Total_customer)
                    Total_customer = dd2;
            }
            myReader002.Close();

            for (int ccc = 1; ccc < Total_customer + 1; ccc++)
            {
                Total_accept_customer = Total_accept_customer + Status_customer[ccc];

                pickup_WD_total = pickup_WD_total + pickup_WD[ccc];
                dropoff_WD_total = dropoff_WD_total + dropoff_WD[ccc];

                if (pickup_WD[ccc] > 0 || dropoff_WD[ccc] > 0)
                {
                    Fare_FT = Fare_FT + g_truej[ccc];
                }
            }

            Total_customer_control[all_control] = Total_customer;
            Acceptance_control[all_control] = Total_accept_customer;

            pickup_WD_total_control[all_control] = pickup_WD_total;
            dropoff_WD_total_control[all_control] = dropoff_WD_total;

            Fare_FT_control[all_control] = Fare_FT;
            //=========================================================================

            objConnection.Close();
        }

        public void Floating_target_results_Excel_Table0921_220421_220422(int all_control, int TYPE)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
   
       
    
            {
                if (TYPE == 1)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                }
                if (TYPE == 9)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                }
                if (TYPE == 3)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                }
                if (TYPE == 5)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                }
                if (TYPE == 6)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                }
                if (TYPE == 4)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }

                if (TYPE == 1201)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1202)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1203)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1204)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1205)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1206)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1207)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1208)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1209)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1210)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1211)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1212)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1213)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1214)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }
            }

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

  
            int[] Status_customer = new int[20000];
            int Total_customer = 0;
            int Total_accept_customer = 0;

            double accept_Q_times_request_ALL_total = 0;
            double accept_Q_times_request_FT_total = 0;
            double accept_Q_times_request_pickup_FT_total = 0;
            double accept_Q_times_request_dropoff_FT_total = 0;
            double Q_times_request_ALL_total = 0;


            double[] pickup_WD = new double[20000];
            double[] dropoff_WD = new double[20000];
            double pickup_WD_total = 0;
            double dropoff_WD_total = 0;
            double Fare_FT = 0;


            string str2 = "select * from customer_status";
            OleDbCommand myCommand002 = new OleDbCommand(str2, objConnection);
            OleDbDataReader myReader002 = myCommand002.ExecuteReader();

            while (myReader002.Read())
            {
                int dd2 = myReader002.GetInt32(1);
                Status_customer[dd2] = myReader002.GetInt16(2);

                pickup_WD[dd2] = myReader002.GetDouble(6);
                dropoff_WD[dd2] = myReader002.GetDouble(14);

                if (dd2 > Total_customer)
                    Total_customer = dd2;
            }
            myReader002.Close();

            for (int ccc = 1; ccc < Total_customer + 1; ccc++)
            {
                Total_accept_customer = Total_accept_customer + Status_customer[ccc];
                accept_Q_times_request_ALL_total = accept_Q_times_request_ALL_total + Status_customer[ccc] * q_truej[ccc];
                Q_times_request_ALL_total = Q_times_request_ALL_total + q_truej[ccc];

                pickup_WD_total = pickup_WD_total + pickup_WD[ccc];
                dropoff_WD_total = dropoff_WD_total + dropoff_WD[ccc];

                if (pickup_WD[ccc] > 0 || dropoff_WD[ccc] > 0)
                {
                    Fare_FT = Fare_FT + g_truej[ccc];

                    accept_Q_times_request_FT_total = accept_Q_times_request_FT_total + Status_customer[ccc] * q_truej[ccc];
                }

                //========================================================================================
                if (pickup_WD[ccc] > 0)
                {
                    accept_Q_times_request_pickup_FT_total = accept_Q_times_request_pickup_FT_total + Status_customer[ccc] * q_truej[ccc];
                }
                if (dropoff_WD[ccc] > 0)
                {
                    accept_Q_times_request_dropoff_FT_total = accept_Q_times_request_dropoff_FT_total + Status_customer[ccc] * q_truej[ccc];
                }
            }

            Total_customer_control[all_control] = Total_customer;
            Acceptance_control[all_control] = Total_accept_customer;

            pickup_WD_total_control[all_control] = pickup_WD_total;
            dropoff_WD_total_control[all_control] = dropoff_WD_total;


            Fare_FT_control[all_control] = Fare_FT;
            //=========================================================================

            accept_Q_times_request_pickup_FT_total_control[all_control] = accept_Q_times_request_pickup_FT_total;
            accept_Q_times_request_dropoff_FT_total_control[all_control] = accept_Q_times_request_dropoff_FT_total;
            accept_Q_times_request_ALL_total_control[all_control] = accept_Q_times_request_ALL_total;
            accept_Q_times_request_FT_total_control[all_control] = accept_Q_times_request_FT_total;
            Q_times_request_ALL_total_control[all_control] = Q_times_request_ALL_total;

            RATE_accept_Q_times_request_pickup_FT_total_control[all_control] = accept_Q_times_request_pickup_FT_total / Q_times_request_ALL_total;
            RATE_accept_Q_times_request_dropoff_FT_total_control[all_control] = accept_Q_times_request_dropoff_FT_total / Q_times_request_ALL_total;
            RATE_accept_Q_times_request_FT_total_control[all_control] = accept_Q_times_request_FT_total / Q_times_request_ALL_total;
            RATE_accept_Q_times_request_ALL_total_control[all_control] = accept_Q_times_request_ALL_total / Q_times_request_ALL_total;

            objConnection.Close();
        }
     
        public void Floating_target_results_Excel_Table1029_220419(int all_control, int TYPE)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
         
  
            {
                if (TYPE == 1)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                }
                if (TYPE == 9)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                }
                if (TYPE == 3)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                }
                if (TYPE == 5)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                }
                if (TYPE == 6)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                }
                if (TYPE == 4)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }

                if (TYPE == 1201)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1202)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1203)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1204)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1205)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1206)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1207)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1208)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1209)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1210)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1211)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1212)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1213)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1214)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }
            }

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();


            string str2 = "select * from customer_coordinates";
            OleDbCommand myCommand002 = new OleDbCommand(str2, objConnection);
            OleDbDataReader myReader002 = myCommand002.ExecuteReader();

            int[,] home_intersection_vehicle_stop = new int[10000, 10000];
            int[,] FT_intersection_vehicle_stop = new int[10000, 10000];
            int totalvehicle = 0;
            int[] totalstop_vehicle = new int[10000];
            double home_vehicle_miles_driving = 0;
            double FT_vehicle_miles_driving = 0;

            while (myReader002.Read())
            {
                int dd2_0 = myReader002.GetInt16(1);
                if (totalvehicle < dd2_0)
                    totalvehicle = dd2_0;
                int dd2_1 = myReader002.GetInt32(2);
                if (totalstop_vehicle[dd2_0] < dd2_1)
                    totalstop_vehicle[dd2_0] = dd2_1;
                home_intersection_vehicle_stop[dd2_0, dd2_1] = myReader002.GetInt16(6);
                FT_intersection_vehicle_stop[dd2_0, dd2_1] = myReader002.GetInt16(9);
            }

            for (int ccc = 1; ccc < totalvehicle + 1; ccc++)
            {

                for (int ddd = 1; ddd < totalstop_vehicle[ccc] + 1; ddd++)
                {
                    int home_a = home_intersection_vehicle_stop[ccc, ddd - 1];
                    int home_b = home_intersection_vehicle_stop[ccc, ddd];
                    home_vehicle_miles_driving = home_vehicle_miles_driving + Intersection_car_d_ij[home_a, home_b];

                    int FT_a = FT_intersection_vehicle_stop[ccc, ddd - 1];
                    int FT_b = FT_intersection_vehicle_stop[ccc, ddd];
                    FT_vehicle_miles_driving = FT_vehicle_miles_driving + Intersection_car_d_ij[FT_a, FT_b];
                }
            }
            home_vehicle_miles_driving_control[all_control] = home_vehicle_miles_driving;
            FT_vehicle_miles_driving_control[all_control] = FT_vehicle_miles_driving;

            objConnection.Close();
        }

        public void Floating_target_results_Excel_Table1_220322_220419(int all_control, int TYPE)
        {
            strConnection = "Provider=Microsoft.ACE.OLEDB.12.0;";
   
       

            {
                if (TYPE == 1)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                }
                if (TYPE == 9)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                }
                if (TYPE == 3)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                }
                if (TYPE == 5)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                }
                if (TYPE == 6)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                }
                if (TYPE == 4)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                    if (all_control == 2)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                    if (all_control == 7)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                    if (all_control == 8)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                    if (all_control == 10)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                    if (all_control == 11)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                    if (all_control == 12)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                    if (all_control == 13)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                    if (all_control == 14)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }

                if (TYPE == 1201)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1201\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1202)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1202\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1203)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1203\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1204)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1204\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1205)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1205\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1206)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1206\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1207)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1207\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1208)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1208\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1209)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1209\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1210)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1210\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1211)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1211\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1212)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1212\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1213)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1213\Output_Instance_file_v2_4.mdb";
                }
                if (TYPE == 1214)
                {
                    if (all_control == 1)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_1.mdb";
                    if (all_control == 9)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_9.mdb";
                    if (all_control == 3)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_3.mdb";
                    if (all_control == 5)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_5.mdb";
                    if (all_control == 6)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_6.mdb";
                    if (all_control == 4)
                        strConnection += @"Data Source=|DataDirectory|\1214\Output_Instance_file_v2_4.mdb";
                }
            }

            OleDbConnection objConnection = new OleDbConnection(strConnection);
            objConnection.Open();

            //=========================================================================
            double Total_profit = 0;
            int Total_vehicle = 0;

            string str1 = "select * from vehicle_profits";
            OleDbCommand myCommand001 = new OleDbCommand(str1, objConnection);
            OleDbDataReader myReader001 = myCommand001.ExecuteReader();

            while (myReader001.Read())
            {
                int dd1 = myReader001.GetInt16(1);
                accum_profit_v[dd1] = myReader001.GetDouble(2); 

                if (dd1 > Total_vehicle)
                    Total_vehicle = dd1;
            }
            myReader001.Close();

            for (int vvv = 1; vvv < Total_vehicle + 1; vvv++)
            {
                Total_profit = Total_profit + accum_profit_v[vvv];
            }
            Profit_control[all_control] = Total_profit;
            //=========================================================================
            int[] Status_customer = new int[20000];
            int[] FT_pickup_customer = new int[20000];
            int[] FT_dropoff_customer = new int[20000];
            int Total_customer = 0;
            int Total_accept_customer = 0;
            int Total_FT_pickup = 0;
            int Total_FT_dropoff = 0;
            int Total_FT_All;

            string str2 = "select * from customer_status";
            OleDbCommand myCommand002 = new OleDbCommand(str2, objConnection);
            OleDbDataReader myReader002 = myCommand002.ExecuteReader();

            while (myReader002.Read())
            {
                int dd2 = myReader002.GetInt32(1);
                Status_customer[dd2] = myReader002.GetInt16(2);
                FT_pickup_customer[dd2] = myReader002.GetInt16(5);
                FT_dropoff_customer[dd2] = myReader002.GetInt16(13);

                if (dd2 > Total_customer)
                    Total_customer = dd2;
            }
            myReader002.Close();

            for (int ccc = 1; ccc < Total_customer + 1; ccc++)
            {
                Total_accept_customer = Total_accept_customer + Status_customer[ccc];
                Total_FT_pickup = Total_FT_pickup + FT_pickup_customer[ccc];
                Total_FT_dropoff = Total_FT_dropoff + FT_dropoff_customer[ccc];
            }
            Total_FT_All = Total_FT_pickup + Total_FT_dropoff;

            Total_customer_control[all_control] = Total_customer;
            Acceptance_control[all_control] = Total_accept_customer;
            Acceptance_rate_control[all_control] = ((double)Total_accept_customer / (double)Total_customer);
            All_Floating_control[all_control] = Total_FT_All;
            FT_Pickup_control[all_control] = Total_FT_pickup;
            FT_Dropoff_control[all_control] = Total_FT_dropoff;
            All_Floating_rate_control[all_control] = ((double)Total_FT_All / ((double)Total_accept_customer * 2));
            FT_Pickup_rate_control[all_control] = ((double)Total_FT_pickup / ((double)Total_accept_customer * 2));
            FT_Dropoff_rate_control[all_control] = ((double)Total_FT_dropoff / ((double)Total_accept_customer * 2));
            //=========================================================================
            int Total_subpath = 0;
            int Total_reposition_subpath = 0;
            int Total_noreposition_subpath = 0;
            int Total_share_subpath = 0;
            int Total_unshare_subpath = 0;
            int[] Share_subpath = new int[100000];

            string str3 = "select * from subpath_profits";
            OleDbCommand myCommand003 = new OleDbCommand(str3, objConnection);
            OleDbDataReader myReader003 = myCommand003.ExecuteReader();

            while (myReader003.Read())
            {
                int dd3 = myReader003.GetInt32(1);
                Share_subpath[dd3] = myReader003.GetInt16(10);

                if (dd3 > Total_subpath)
                    Total_subpath = dd3;
            }
            myReader003.Close();

            for (int ppp = 1; ppp < Total_subpath + 1; ppp++)
            {
                if (Share_subpath[ppp] == 0)
                {
                    Total_reposition_subpath = Total_reposition_subpath + 1;
                }
                else
                {
                    Total_noreposition_subpath = Total_noreposition_subpath + 1;
                    if (Share_subpath[ppp] == 1)
                        Total_unshare_subpath = Total_unshare_subpath + 1;
                    if (Share_subpath[ppp] > 1)
                        Total_share_subpath = Total_share_subpath + 1;
                }
            }

            Total_unreposition_path_control[all_control] = Total_noreposition_subpath;
            Shared_path_control[all_control] = Total_share_subpath;
            Shared_path_rate_control[all_control] = ((double)Total_share_subpath / (double)Total_noreposition_subpath);

            objConnection.Close();

        }


        public void close_writer()
        {
            writer.Close();
        }

        static void Main(string[] args)
        {
            DateTime begin_total;
            DateTime end_total;
            TimeSpan interval_total;
            begin_total = DateTime.Now;

            Program prg = new Program();
            prg.Rolling_Horizon_overall0202_0212_flexibledropoff_matching_clearversion_reposition_changej(50, W, 4);
            prg.close_writer();
            Console.Read();

            //====================================================================================================== Table 3
            prg = new Program();
            prg.Floating_target_results_Excel_Table0524_ALL_220419(4);
            prg.close_writer();

            prg = new Program();
            prg.Floating_target_results_Excel_Table1_ALL_220322_220419(4);
            prg.close_writer();
    
            prg = new Program();
            prg.Floating_target_results_Excel_Table1029_ALL0_220419(4);
            prg.close_writer();
            Console.Read();
            //==============================================================================================      

            end_total = DateTime.Now;
            interval_total = end_total.Subtract(begin_total);
            Console.WriteLine("T_0 = " + T_0 + " customer = " + n + ";" + " Q = " + Q + "; apha = " + apha + ";");
            Console.WriteLine("*******************************************************************");
            Console.WriteLine("DONE! The Total Time is " + Convert.ToString(interval_total.Hours) + ":" + Convert.ToString(interval_total.Minutes) + ":" + Convert.ToString(interval_total.Seconds));
  
        }
    }
}



