using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace Honeywell.UOP.CPP2CS
{
    public struct UserDataStruct
    {
        public int n;
        public double x;
    }

    public enum TemperatureProfileEnum
    {
        ISOTHERMAL = 0,
        ADIABATIC = 1,
        PROFILE = 2
    };

    class Program
    {
        static void Main(string[] args)
        {
            int n = IntPtr.Size;
            int size = 2;
            int[] a = { 1, 2 };
            int[] b = { 3, 4 };
            PInvoker.TestMPI(size, a, b);

            //PInvoker.test3(size, a, b);

            TemperatureProfileEnum reactorTemperatureProfile = TemperatureProfileEnum.PROFILE;
            double[] yin = { 1.0, 2.0, 3.0 };
            double[] yout = { 0.0, 0.0, 0.0 };
            double[] temperatureProfile = null;
            PInvoker.RunUnivModel(reactorTemperatureProfile, 2, yin, yout, 0, temperatureProfile);

        }
    }

}
