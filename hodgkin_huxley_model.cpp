/*
    Hodgkin-Huxley Model
        A continuous-time dynamical system model of biophysically based neuron to explain the initiation and propagation of action potentials
        through ion channels of a membrane.

        It consists of a set of non-linear differential equations
            Membrane Capacitance
                Cm dV / dt = gNa m^3h(E_Na - V) + gK n^4(E_K - V) + gL(E_L - V)
            
            Potassium Gating Variable: n
                dm / dt = alpha_m(V)(1-m) - beta_m(V)m

            Sodium Activation Gating Varible: m
                dh / dt = alpha_h(V)(1-m) - beta_h(V)h

            Sodium Inactvation Gating Varible: h
                dn / dt = alpha_n(V)(1-n) - beta_n(V)n      

    References
        The Core Equation of Neuroscience: https://www.youtube.com/watch?v=zOmhHE2xctw
        Computational Neuroscience:        https://mrgreene09.github.io/computational-neuroscience-textbook/Ch4.html
        Neuronal Dynamics:                 https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
        Hodgkin-Huxley Tutorial:           http://nelson.beckman.illinois.edu/courses/physl317/part1/Lec3_HHsection.pdf

    Compile
        g++ ./hodgkin_huxley_model.cpp -o ./hodgkin_huxley_model.exe
        ./hodgkin_huxley_model.exe
*/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>

// Parameters
// Membrane capacitance: uF/cm^2
constexpr float Cm = 1.0f;    

// Maximal conductance: mS/cm^2
constexpr float gK  = 35.0f; // maximal K conductance
constexpr float gNa = 40.0f; // maximal Na conductance
constexpr float gL  = 0.3f;  // maximal leak conductance

// Nernst Equilibrium Potential: mV
constexpr float EK  = -77.0f; // K conductance reversal potential
constexpr float ENa = 55.0f;  // Na conductance reversal potential
constexpr float EL  = -65.0f; // Leak conductance reversal potential

// Membrane internal current
float potassium_channel(float V, float n) { return gK * std::pow(n, 4) * (V - EK); }
float sodium_channel(float V, float m, float h) { return gNa * std::pow(m, 3) * h * (V - ENa); }
float leak_channel(float V) { return gL * (V - EL); }
float I_int(float V, float n, float m, float h) { return potassium_channel(V, n) + sodium_channel(V, m, h) + leak_channel(V); }

// Gating variables
float alpha_n(float V) { return 0.01f  * (V + 55.0f) / (1.0f - std::exp(-0.1f * (V + 55.0f))); }
float beta_n(float V)  { return 0.125f * std::exp(-0.01125f * (V + 55.0f)); }
float alpha_m(float V) { return 0.1f   * (V + 40.0f) / (1.0f - std::exp(-0.1f * (V + 40.0f)));}
float beta_m(float V)  { return 4.0f   * std::exp(-0.05556f * (V + 65.0f)); }
float alpha_h(float V) { return 0.07f  * std::exp(-0.05f * (V + 65.0f)); }
float beta_h(float V)  { return 1.0f   / (1.0f + std::exp(-0.1f * (V + 35.0f))); }

// Membrane external current
float I_ext(float t) { return (t > 10.0f && t < 40.0f) ? 10.0f : 0.0f; }

int main()
{
    // Initial conditions
    float T  = 50.0f;
    float dt = 0.01f;

    float V  = -65.0f;  // membrane voltage
    float n  = 0.3177f; // potassium gate
    float m  = 0.0529f; // sodium activation gate
    float h  = 0.5961f; // sodium inactivation gate

    std::ofstream out("hh_data.csv");
    out << "t\tV\tn\tm\th\n";

    for (float t = 0.0f; t < T; t += dt)
    {
        // Compute the gate permissive or non-permissive state
        float dn = (alpha_n(V) * (1.0f - n) - beta_n(V) * n);
        float dm = (alpha_m(V) * (1.0f - m) - beta_m(V) * m);
        float dh = (alpha_h(V) * (1.0f - h) - beta_h(V) * h);

        n += dn * dt;
        m += dm * dt;
        h += dh * dt;

        // Compute voltage difference
        float dV = (I_ext(t) - I_int(V, n, m, h)) / Cm;

        V += dV * dt;

        out << t << "\t" << V << "\t" << n << "\t" << m << "\t" << h << "\n";
    }

    out.close();

    return 0;
}