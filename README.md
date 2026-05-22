# Capstone
Code from my ECE Capstone project at Carnegie Mellon University

18-500 Spring 2025

This code controls a microprocessor for a digital effects pedal which uses customizable sequence patterns, pitch modulation, tempo, and delay attenuation.

Code was written in C++ for the Daisy Seed microprocessor, in collaboration with teammates Nicholas Walker and Chaitanya Irkar. Chaitanya worked primarily on the physical circuitry, Nicholas worked primarily on the digital-analog interfacing, and I worked primarily on the signal chain algorithms, including for pitch modulation via phase vocoder.

Schedule and hardware constraints meant that we had to be clever with some of our software. The microcontroller wasn't cooperating with importing libraries for the discrete Fourier transform, so I reprogrammed it myself with the course notes I had taken from 18-491.

Functions that I worked on:
- (...)
