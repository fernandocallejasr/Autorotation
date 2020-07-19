import pandas as pd
import numpy as np

class ElemPala:
    """Clase que te permite obtener los parámetros diferenciales (Levantamiento, Arrastre, Torque)
       de cada elemento de una pala de un autogiro o helicóptero, así como el torque total de cada 
       elemento en las dos palas. También se puede obtener el número de revoluciones necesario para
       entrar en autorotación cuando el torque total es cero
    """
    
    def __init__(self, Data, omega, dr, hdens):
        self.omega = omega
        
        self.Data = Data
        self.W = Data['W']
        self.D = Data['D'] 
        self.cuerda = Data['cuerda']; 
        self.a = Data['a'] 
        self.Cdo = Data['Cdo']
        self.long_efec = Data['long_efec']
        self.alpha_rotor = Data['alpha_rotor']
        self.Voo = Data['Voo']
        self.theta = Data['theta']
        
        self.dr = dr
        self.num_elem = int(self.long_efec/self.dr)
        self.r_0 = self.D/2 - self.long_efec
        
        self.hdens = hdens
        self.dens   = 1.225*(1-self.hdens*(0.00357/518.4))**(1/0.235)
        self.dens_ratio = self.dens/1.225
        
        self.T = self.W/np.cos(np.deg2rad(self.alpha_rotor))  
    
        
        self.list_r = [0 for x in range(self.num_elem)]
        self.list_alpha_fb = [0 for x in range(self.num_elem)] # alphas de los elementos de la pala que avanza
        self.list_alpha_bb = [0 for x in range(self.num_elem)] # alphas de los elementos de la para que retrocede
        self.list_VR_fb = [0 for x in range(self.num_elem)]; self.list_VR_bb = [0 for x in range(self.num_elem)]
        self.list_L_fb = [0 for x in range(self.num_elem)]; self.list_L_bb = [0 for x in range(self.num_elem)] 
        self.list_D_fb = [0 for x in range(self.num_elem)]; self.list_D_bb = [0 for x in range(self.num_elem)]
        self.list_Q_fb = [0 for x in range(self.num_elem)]; self.list_Q_bb = [0 for x in range(self.num_elem)]
        
        self.DF = pd.DataFrame({'Omega [rpm]': [0], 'Q total': [0]})
        self.Q = []
        
    def calc_radios(self):
        """
        Método para obtener el radio correspondiente a cada elemento
        Args: None
        Output: None
        """
        self.list_r = [round(r,2) for r in np.arange(self.r_0 + self.dr, self.D/2, self.dr)]
        
    def calc_VR(self):
        """Método para obtener VR en cada elemento para la pala que avanza y
           la que retrocede
           Args: None
           Output: None
        """
        i = 0
        alpha = np.deg2rad(self.alpha_rotor)
        for r in self.list_r:
            self.list_VR_fb[i] = ((self.Voo*np.sin(alpha))**2 + (self.omega*r*(1/60)*(2*np.pi) + self.Voo*np.cos(alpha))**2)**0.5
            self.list_VR_bb[i] = ((self.Voo*np.sin(alpha))**2 + (self.omega*r*(1/60)*(2*np.pi) - self.Voo*np.cos(alpha))**2)**0.5
            i += 1
        
    def calc_alpha(self):
        """Método que nos permite obtener el ángulo de ataque de cada elemento en radianes
           Args: None
           Output: None
        """
        i = 0
        for r in self.list_r:
            self.list_alpha_fb[i] = np.arcsin((self.Voo * np.sin(np.deg2rad(self.alpha_rotor)))/self.list_VR_fb[i]) + self.theta*(np.pi/180)
            self.list_alpha_bb[i] = np.arcsin((self.Voo * np.sin(np.deg2rad(self.alpha_rotor)))/self.list_VR_bb[i]) + self.theta*(np.pi/180)
            i += 1
    
    def calc_delta_L(self):
        """Método para calcular el delta de Levantamiento en cada elemento
        """
        i = 0
        for VR in self.list_VR_fb:
            self.list_L_fb[i] = 0.5*self.dens*(VR**2)*self.cuerda*self.a*self.list_alpha_fb[i]*self.dr
            i += 1
            
        i = 0
        for VR in self.list_VR_bb:
            self.list_L_bb[i] = 0.5*self.dens*(VR**2)*self.cuerda*self.a*self.list_alpha_bb[i]*self.dr
            i += 1
            
    def calc_delta_D(self):
        """Método para calcular el delta de Levantamiento en cada elemento
        """
        i = 0
        for VR in self.list_VR_fb:
            self.list_D_fb[i] = 0.5*self.dens*(VR**2)*self.cuerda*self.Cdo*self.dr
            i += 1
            
        i = 0
        for VR in self.list_VR_bb:
            self.list_D_bb[i] = 0.5*self.dens*(VR**2)*self.cuerda*self.Cdo*self.dr
            i += 1
        
    
    def calc_torque(self):
        """Método para calcular el torque total de ambas palas
        """
        for i in range(len(self.list_VR_fb)):
            self.list_Q_fb[i] = (self.list_D_fb[i] - self.list_alpha_fb[i]*self.list_L_fb[i])*self.list_r[i]*self.dr
            self.list_Q_bb[i] = (self.list_D_bb[i] - self.list_alpha_bb[i]*self.list_L_bb[i])*self.list_r[i]*self.dr
        
        self.Q = sum(self.list_Q_fb) + sum(self.list_Q_bb)
    
    def calc_All(self):
        """Método que calcula todos los parámetros necesarios para obtener el torque
        """
        self.__init__(self.Data, self.omega, self.dr, self.hdens)
        self.calc_radios()
        self.calc_VR()
        self.calc_alpha()
        self.calc_delta_L()
        self.calc_delta_D()
        self.calc_torque()
    
    def get_DataFrame_All(self):
        """Método que calcula todos los parámetros de todos los elementos de la pala y los pone en un Data Frame de pandas
        """
        self.calc_All()
        
        DF_total = pd.DataFrame({'r [m]': self.list_r,
                                'VR fb [m/s]': self.list_VR_fb,
                                'VR bb [m/s]': self.list_VR_bb,
                                'alpha fb [°]': np.degrees(self.list_alpha_fb),
                                'alpha bb [°]': np.degrees(self.list_alpha_bb),
                                'dL fb [N]': self.list_L_fb,
                                'dL bb [N]': self.list_L_bb,
                                'dD fb [N]': self.list_D_fb,
                                'dD bb [N]': self.list_D_bb,
                                'dQ fb [N]': self.list_Q_fb,
                                'dQ bb [N]': self.list_Q_bb 
                               })
        
        return DF_total
    
    def calc_list_torque(self, list_omega):
        """Método que calcula el torque para una lista de distintas rpm
           Inputs:
                    list_omega: lista. Lista de diferentes rpm para calcular su torque
           Outputs:
                   Data Frame de rpm y su respectivo torque
        """
        list_Q = [0 for x in range(len(list_omega))]
        
        i = 0
        for omega in list_omega:
            self.omega = omega
            self.calc_All()
            list_Q[i] = self.Q
            i += 1
            
        DF_torque = pd.DataFrame({'RPM': omega, 'Torque [Nm]': list_Q})
        return DF_torque
            
    def calc_rpm(self, rpm_inicial):
        """Método que varía las rpm hasta encontrar las que dan el estado de autorotación
           Inputs:
               rpm_inicial (float): RPM para comenzar a iterar 
           Outputs:
               RPM para autorotación
        """
        Error = 1
        rpm = rpm_inicial
        self.omega = rpm
        
        while Error > 0.01:
            self.omega = rpm
            self.calc_All()
            Error = abs(self.Q)
            rpm += 0.001
        
        return f'RPM para autorotación: {round(rpm, 4)}'

def main():
    AG1_Data = {'W': 450*9.81,      #[N]
               'D' : 5,           #[m]
               'cuerda': 0.24,      #[m]
               'a': 4.8,            #[/rad]
               'Cdo' : 0.015,       #[1]
               'long_efec': 3.68,   #[m]
               'alpha_rotor': 30,   
               'Voo_Vh': 2, 
               'theta': 0
                }
    Final = ElemPala(AG1_Data, 546, 0.08, 0)
    Final.calc_rpm(530)

import unittest

AG1_Data = {'W': 450*9.81,
           'D' : 8.4,
           'cuerda': 0.24,
           'a': 5.6,
           'Cdo' : 0.007, 
           'long_efec': 3.68,
           'alpha_rotor': 30,
           'Voo_Vh': 2, 
           'theta': 0
            }

class Test_bemClass(unittest.TestCase):
    def setUp(self):
        self.bem = ElemPala(AG1_Data, 500, 0.08, 0)
        self.bem.calc_All()
        
    def test_initializaion(self):
        #self.assertEqual(self.clothing.color, 'orange', 'color should be orange')
        self.assertEqual(round(self.bem.T,2), 5097.43, 'T = 5097.43 N')
        self.assertEqual(self.bem.list_r[4] , 0.92, 'Radio 5° Elem = 0.92 m')
        self.assertEqual(round(self.bem.Voo,2), 12.25, 'Voo = 12.25 m/s')
        self.assertEqual(round(self.bem.list_alpha_fb[4], 4) , 0.1039, 'Alpha 5° Elem fb = 0.1038 rad')
        self.assertEqual(round(self.bem.list_alpha_bb[4], 4) , 0.1617, 'Alpha 5° Elem bb = 0.1617 rad')
        self.assertEqual(round(self.bem.list_L_fb[4], 2) , 23.64, 'L 5° Elem fb = 23.64 N')
        self.assertEqual(round(self.bem.list_L_bb[4], 2) , 15.02, 'L 5° Elem bb = 15.02 rad')
        self.assertEqual(round(self.bem.list_D_fb[4], 2) , 0.28, 'D 5° Elem fb = 0.28 N')
        self.assertEqual(round(self.bem.list_D_bb[4], 2) , 0.12, 'D 5° Elem bb = 0.12 N')
        self.assertEqual(round(self.bem.list_Q_fb[4], 2) , -0.16, 'Q 5° Elem fb = -0.16 Nm')
        self.assertEqual(round(self.bem.list_Q_bb[4], 2) , -0.17, 'Q 5° Elem bb = -0.17 Nm')
        
tests = Test_bemClass()
tests_loaded = unittest.TestLoader().loadTestsFromModule(tests)
unittest.TextTestRunner().run(tests_loaded)

if __name__ == '__main__':
  main()