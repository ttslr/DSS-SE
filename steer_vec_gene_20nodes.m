
clear;clc 
% steering vector generation
 Num = 20;
 
 s = [6,4];
%  s = [6.5,5];
%  s = [4.5,3];
%  s = [2,5];
%  s = [3,2];
 
 
 
 mic(1,:) = [2,2];
 mic(2,:) = [1.2,3.7];
 mic(3,:) = [1.7,5.4];
 mic(4,:) = [2.6,7.1];
 mic(5,:) = [3.6,1.4];
 mic(6,:) = [3.2,3.3];
 mic(7,:) = [3.8,4.1];
 mic(8,:) = [3.5,6.8];
 mic(9,:) = [5.1,1.9];
 mic(10,:) = [5.5,3.5];
 mic(11,:) = [5.0,5.1];
 mic(12,:) = [5.6,6.6];
 mic(13,:) = [6.3,1.5];
 mic(14,:) = [6.6,3.4];
 mic(15,:) = [6.4,5.2];
 mic(16,:) = [6.5,7.3];
 mic(17,:) = [8.2,1.6];
 mic(18,:) = [8.1,4.2];
 mic(19,:) = [8.1,4.8];
 mic(20,:) = [7.6,6.8];
 
 for i = 1:Num
     avec_a_20nodes(i,1) = 1/norm(s-mic(i,:),2);
     avec_tao_20nodes(i,1) = norm(s-mic(i,:),2)/343;     
 end
 
 save('avec_a_20nodes','avec_a_20nodes');
 save('avec_tao_20nodes','avec_tao_20nodes');
 
 
  for i = 1:Num
      for j = 1:Num
          Dmatrix(i,j) = norm( (mic(i,:)-mic(j,:)) ,2 );
      end  
  end
 
   save('Dmatrix','Dmatrix');