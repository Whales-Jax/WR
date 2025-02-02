#pragma once
/*
class parameter
{
private:
	//refrigerant
	struct data{
		double T_in_r;			//refrigerant temp in C
		double T_out_r;			//refrigerant temp out C
		double P_in_r;			//refrigerant pressure in kPa
	};

	data **param;
public:
	const int length;
	parameter(int length);
	parameter();
	~parameter();
};

parameter::parameter(int length) : length(length) {
		param = new data*[length];
		for (int i = 0; i < length; i++) {
				param[i] = new data[length];
		}
}
*/
class parameter
{
public:

	//refrigerant
	struct refr{
		//===Refrigerant Local Properties at previous time step===//
		double T_R_P_dt;												//refrigerant local temperature [degC]
		double P_R_P_dt;												//refrigerant local pressure [kPa]
		double h_R_P_dt;												//refrigerant local specific enthalpy [kJ/kg]
		double u_R_P_dt;												//refrigerant local specific internal energy [kJ/kg]
		double M_R_P_dt;

		//===Refrigerant Boundary Properties at previous time step===//
		double G_R_B_dt;												//refrigerant boundary mass flowrate [kg/s]
		double P_R_B_dt;												//refrigerant boundary pressure [kPa]
		double h_R_B_dt;												//refrigerant boundary specific enthalpy [kJ/kg]
		double T_R_B_dt;												//refrigerant boundary temperature [degC]
	
		//===Refrigerant Local Properties current time step===//
		double T_R_P;													//refrigerant local temperature [degC]
		double P_R_P;													//refrigerant local pressure [kPa]
		double h_R_P;													//refrigerant local specific enthalpy [kJ/kg]
		double u_R_P;													//refrigerant local specific internal energy [kJ/kg]
		double M_R_P;
		double rho_R_P;													//refrigerant local density [kg/m^3]
		double myu_R_P;													//refrigerant local viscosity [Pa-s]
		double k_R_P;													//refrigerant local thermal conductivity [W/m-k]
		double cp_R_P;													//refrigerant local specific heat [J/kg-K]
		double x_R_P;													//refrigerant local quality
		double VF_R_P;													//refrigerant local void fraction
		double myu_Rl_P;
		double k_Rl_P;
		double rho_Rl_P;
		double cp_Rl_P;
		double h_Rl_P;
		double myu_Rv_P;
		double k_Rv_P;
		double rho_Rv_P;
		double cp_Rv_P;
		double h_Rv_P;
		double st_Rlv_P;

		double Re_R_P;													//
		double Pr_R_P;													//

		//===Refrigerant Boundary Properties current time step===//
		double G_R_B;													//refrigerant boundary mass flowrate [kg/s]
		double P_R_B;													//refrigerant boundary pressure [kPa]
		double h_R_B;													//refrigerant boundary specific enthalpy [kJ/kg]
		double T_R_B;													//refrigerant boundary temperature [degC]

		double rho_R_B;													//refrigerant boundary density [kg/m^3]
		double myu_R_B;													//refrigerant boundary viscosity [Pa-s]
		double k_R_B;													//refrigerant boundary thermal conductivity [W/m-k]
		double cp_R_B;													//refrigerant boundary specific heat [J/kg-K]
		double x_R_B;													//refrigerant boundary quality

		double Re_R_B;													//
		double Pr_R_B;	
	
	};

	struct air {
		//===Air Local Properties at previous time step===//
		double P_A_P_dt;												//air local ipressure [kPa]
		double h_A_P_dt;												//air local specific enthalpy [kJ/kg]
		double T_A_P_dt;												//air local temperature [degC]
		double u_A_P_dt;												//air local specific internal energy [kJ/kg]
		double v_A_P_dt;												//air local velocity [m/s]

		//===Air Boundary Properties at previous time step===//
		double G_A_B_dt;												//air boundary mass flowrate [kg/s]
		double h_A_B_dt;												//air boundary specific enthalpy [kJ/kg]
		double P_A_B_dt;												//air boundary pressure [kPa]
		double T_A_B_dt;												//air boundary temperature [degC]
		double v_A_B_dt;
	


		//===Air Local Properties current time step===//
		double P_A_P;													//air local ipressure [kPa]
		double h_A_P;													//air local specific enthalpy [kJ/kg]
		double T_A_P;													//air local temperature [degC]
		double u_A_P;													//air local specific internal energy [kJ/kg]
		double v_A_P;													//air local velocity [m/s]
		double rho_A_P;													//air local density [kg/m^3]
		double myu_A_P;													//air local viscosity [Pa-s]
		double k_A_P;													//air local thermal conductivity [W/m-k]
		double cp_A_P;													//air local specific heat [J/kg-K]
		double x_A_P;													//air local absolute humidity

		double Re_A_P;													//
		double Pr_A_P;													//
		double Remax_A_P;

		//===Air Boundary Properties current time step===//
		double G_A_B;													//air boundary mass flowrate [kg/s]
		double h_A_B;													//air boundary specific enthalpy [kJ/kg]
		double P_A_B;													//air boundary pressure [kPa]
		double T_A_B;													//air boundary temperature [degC]
		double v_A_B;													//air boundary velocity [m/s]
		double u_A_B;													//air boundary specific internal energy [kJ/kg]
		double rho_A_B;													//air boundary density [kg/m^3]
		double myu_A_B;													//air boundary viscosity [Pa-s]
		double k_A_B;													//air boundary thermal conductivity [W/m-k]
		double cp_A_B;													//air boundary specific heat [J/kg-K]
		double x_A_B;													//air boundary absolute humidity

		double Re_A_B;													//
		double Pr_A_B;													//
	};
	struct prop {
		//===Transfer properties at previous time step===//
		double q_P_dt;													//local heat transfer rate [kW]
		double K_P_dt;													//local overall heat transfer coefficeint [kW/m^2-K]

		//===Transfer properties at current time step===//
		double q_P;														//local heat transfer rate [kW]
		double K_P;														//local overall heat transfer coefficeint [kW/m^2-K]


	};
	//refr **param;


	parameter();
	~parameter();
};