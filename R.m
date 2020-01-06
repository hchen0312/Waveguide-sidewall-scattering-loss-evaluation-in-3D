function [ R ] = R( sigma,Lc,x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
R=2*sigma*sigma*Lc/(1+Lc*Lc*x*x);
end

