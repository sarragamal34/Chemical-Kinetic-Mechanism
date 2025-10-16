function [k]=K_Arrhenius(Table_Data,T)
R=0.00198720425864083;
k=zeros(size(Table_Data.A));
for i=1:1:size(k)
A  = Table_Data.A(i) ;               % column A, row i [web:55]
n  = Table_Data.m(i);                % column m, row i [web:55]
Ea = Table_Data.Ea(i);               % column Ea, row i [web:55]
k(i) = A .* T.^n .* exp(-Ea./(R.*T)); % numeric operations [web:55]
end
end