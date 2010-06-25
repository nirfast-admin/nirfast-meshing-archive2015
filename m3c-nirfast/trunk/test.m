for i=1:5
    a=ones(255,255,'uint8')*200;
    a(1:2,:)=0;
    a(:,1:2)=0;
    a(end-2:end,:)=0;
    a(:,end-2:end)=0;
    a(50:80,50:80)=10;
    imwrite(a,['foo' num2str(i) '.bmp'],'bmp');
end

MMC('foo', 2, 4, 4, 'myoutput')