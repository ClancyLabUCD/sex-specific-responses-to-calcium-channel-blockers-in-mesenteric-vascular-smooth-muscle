% x = mean(V)
% a = mean (Ca_in)
% b = mean (Na_in)
% c = mean (K_in)
% d = mean(I_KvALL)
% e = mean (I_Kv15)
% f = mean (I_Kv21)
% g = mean (I_CaL) 


data1 = load ('all_ODEs.txt');
id0 = find( data1(:,1) > 400000 );
zDFQsAvg = [ 
    mean( data1(id0,2) )
    mean( data1(id0,10) ) 
%     mean( data1(id0,8) )
%     mean( data1(id0,10) )
    ]';

% id0 = find( data1(:,1) > 400000 );
% zDFQsAvg = [ 
%     mean( data1(id0,2) )
%     mean( data1(id0,7) ) 
%     mean( data1(id0,8) )
%     mean( data1(id0,10) )
%     ]';
% 
% 
% id1 = find( data2(:,1) > 400000 );
% zCurrentAvg = [ 
%     mean( data2(id1,2) )
%     mean( data2(id1,3) ) 
%     mean( data2(id1,4) )
%     mean( data2(id1,5) )
%     mean( data2(id1,6) )
%     mean( data2(id1,7) ) 
%     mean( data2(id1,8) )
%     mean( data2(id1,9) )
%     mean( data2(id1,10) )
%     mean( data2(id1,11) ) 
%     mean( data2(id1,12) )
%     mean( data2(id1,13) )
%     mean( data2(id1,14) ) 
%     mean( data2(id1,15) )
%     mean( data2(id1,16) )
%     ]';
