clc
clear
q_return=[];
hospital_wait_all=[];
queue=[];
queue_1=[];
queue_all=[];
price=zeros(1,3);
time=zeros(1,3);
H=['A' 'B' 'C'];
u1=1.1;
A=[1 2 3];A2=[1 0 2];A3=[1 2 0];A23=[1 0 0];
B=[2 1 3];B2=[1 0 2];B3=[2 1 0];B23=[1 0 0];
C=[1 3 2];C2=[1 0 2];C3=[1 2 0];C23=[1 0 0];
costA=0;costB=0;costC=0;timeA=0;timeB=0;timeC=0;
kk=[1 2 3];
hospital_pref=[3 2 1; 1 3 2; 2 1 3];
t=.001;tt=0;ttt=0;
hospital_wait=zeros(3,3,3);
patient_visitor=zeros(3,3,3);
Patient_pref=zeros(3,3,3);
time_stamp=[];
time_finish=[];
k_wait=[6 6];
TED=400;
for T=0:1:TED
    for u=8*T:t:8*(T+1)
        lambda=[0.2*t*u1 0.1*t*u1 0.05*t*u1];
        r=poissrnd(lambda);
        if r(1)==1
            queue_all=[queue_all 'A'];
            queue=[queue,[double('A');u]];
        end
        if r(2)==1
            queue_all=[queue_all 'B'];
            queue=[queue,[double('B');u]];
        end
        if  r(3)==1
            queue_all=[queue_all 'C'];
            queue=[queue,[double('C');u]];
        end
    end
    if length(queue(1,:))==3
        n=3;
        n1=3;
    elseif length(queue(1,:))==2
        n=2;
        n1=2;
    elseif length(queue(1,:))==1
        n=1;
        n1=1;
    elseif isempty(queue)
        n=0;
        n1=0;
    else
        n=3;
        n1=length(queue(1,:));
    end
    queue_1(1,:)=queue(1,n+1:n1);
    queue_1(2,:)=queue(2,n+1:n1);
    Patient_pref=zeros(3,3,3);
    if ~isempty(queue) && queue(1,1)==double('A')
        time_stamp=[time_stamp queue(2,1)];
        if kk(2)==0 && kk(3)~=0
            Patient_pref(1,1,:)=A2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(1,1,:)=A3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(1,1,:)=A23;
        else
            Patient_pref(1,1,:)=A;
        end
    elseif ~isempty(queue) && queue(1,1)==double('B')
        if kk(2)==0 && kk(3)~=0
            if k_wait(1)==0
                Patient_pref(2,1,:)=B2;
            else
                k_wait(1)=k_wait(1)-1;
                queue_1(1,:)=[double('B') queue_1(1,2:end)];
                queue_1(2,:)=[queue(2,1) queue_1(2,2:end)];
                queue(1)=[];
                queue(2)=[];
            end
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(2,1,:)=B3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(2,1,:)=B23;
        else
            Patient_pref(2,1,:)=B;
        end
    elseif ~isempty(queue) && queue(1,1)==double('C')
        if kk(2)==0 && kk(3)~=0
            Patient_pref(3,1,:)=C2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(3,1,:)=C3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(3,1,:)=C23;
        else
            Patient_pref(3,1,:)=C;
        end
    end
    if  length(queue(1,:))>1 && queue(1,2)==65
        if kk(2)==0 && kk(3)~=0
            Patient_pref(1,2,:)=A2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(1,2,:)=A3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(1,2,:)=A23;
        else
            Patient_pref(1,2,:)=A;
        end
    elseif length(queue(1,:))>1 && queue(1,2)==double('B')
        if kk(2)==0 && kk(3)~=0
            if k_wait(2)==0
                Patient_pref(2,2,:)=B2;
            else
                k_wait(2)=k_wait(2)-1;
                queue_1(1,:)=[double('B') queue_1(1,2:end)];
                queue_1(2,:)=[queue(2,2) queue_1(2,2:end)];
                queue(3)=[];
                queue(4)=[];
            end
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(2,2,:)=B3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(2,2,:)=B23;
        else
            Patient_pref(2,2,:)=B;
        end
    elseif length(queue(1,:))>1 && queue(1,2)==double('C')
        if kk(2)==0 && kk(3)~=0
            Patient_pref(3,2,:)=C2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(3,2,:)=C3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(3,2,:)=C23;
        else
            Patient_pref(3,2,:)=C;
        end
    end
    if length(queue(1,:))>2 && queue(1,3)==double('A')
        if kk(2)==0 && kk(3)~=0
            Patient_pref(1,3,:)=A2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(1,3,:)=A3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(1,3,:)=A23;
        else
            Patient_pref(1,3,:)=A;
        end
    elseif length(queue(1,:))>2 && queue(1,3)==double('B')
        if kk(2)==0 && kk(3)~=0
            Patient_pref(2,3,:)=B2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(2,3,:)=B3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(2,3,:)=B23;
        else
            Patient_pref(2,3,:)=B;
        end
    elseif length(queue(1,:))>2 && queue(1,3)==double('C')
        if kk(2)==0 && kk(3)~=0
            Patient_pref(3,3,:)=C2;
        elseif kk(2)~=0 && kk(3)==0
            Patient_pref(3,3,:)=C3;
        elseif kk(2)==0 && kk(3)==0
            Patient_pref(3,3,:)=C23;
        else
            Patient_pref(3,3,:)=C;
        end
    end
    for s=1:27
        if Patient_pref(s)==1
            patient_visitor(s)=1;
        end
    end
    kk1=kk;
    kk(kk==0)=[];
    kkk=kk;
    kk=kk1;
    for ll=1:21
        for k=kkk
            f=hospital_pref(k,:);
            if find(patient_visitor(f(1),:,k)==1,1)
                hospital_wait(f(1),find(patient_visitor(f(1),:,k)==1,1),k)=1;
                patient_visitor(f(1),find(patient_visitor(f(1),:,k)==1,1),k)=0;
            elseif find(patient_visitor(f(2),:,k)==1,1)
                hospital_wait(f(2),find(patient_visitor(f(2),:,k)==1,1),k)=1;
                patient_visitor(f(2),find(patient_visitor(f(2),:,k)==1,1),k)=0;
            elseif find(patient_visitor(f(3),:,k)==1,1)
                hospital_wait(f(3),find(patient_visitor(f(3),:,k)==1,1),k)=1;
                patient_visitor(f(3),find(patient_visitor(f(3),:,k)==1,1),k)=0;
            end
        end
        for l=1:27
            if length(kkk)==2
                if patient_visitor(l)==1 && Patient_pref(l)==1
                    [ X, Y, Z]=ind2sub([3 3 3],l);
                    Z1=find(Patient_pref(X,Y,:)==2);
                    hospital_wait(X,Y,Z1)=1;
                    patient_visitor(X,Y,Z)=0;
                elseif patient_visitor(l)==1 && Patient_pref(l)==2
                    [ X, Y, Z]=ind2sub([3 3 3],l);
                    patient_visitor(X,Y,Z)=0;
                    for i=1:3
                        for j=1:3
                            for d=1:3
                                if X==j && Y==i && Z==d
                                    q_return=[q_return ,[double(H(j));queue(2,i)]];
                                    queue(:,i)=[];
                                end
                            end
                        end
                    end
                end
            elseif length(kkk)==1
                if patient_visitor(l)==1 && Patient_pref(l)==1
                    [ X, Y, Z]=ind2sub([3 3 3],l);
                    patient_visitor(X,Y,Z)=0;
                    for i=1:3
                        for j=1:3
                            for d=1:3
                                if X==j && Y==i && Z==d
                                    q_return=[q_return,[double(H(j)); queue(2,i)]];
                                    queue(:,i)=[];
                                end
                            end
                        end
                    end
                end
            else
                if patient_visitor(l)==1 && Patient_pref(l)==1
                    [ X, Y, Z]=ind2sub([3 3 3],l);
                    Z1=find(Patient_pref(X,Y,:)==2);
                    hospital_wait(X,Y,Z1)=1;
                    patient_visitor(X,Y,Z)=0;
                elseif patient_visitor(l)==1 && Patient_pref(l)==2
                    [ X, Y, Z]=ind2sub([3 3 3],l);
                    Z1=find(Patient_pref(X,Y,:)==3);
                    hospital_wait(X,Y,Z1)=1;
                    patient_visitor(X,Y,Z)=0;
                elseif patient_visitor(l)==1 && Patient_pref(l)==3
                    [ X, Y, Z]=ind2sub([3 3 3],l);
                    patient_visitor(X,Y,Z)=0;
                    for i=1:3
                        for j=1:3
                            for d=1:3
                                if X==j && Y==i && Z==d
                                    q_return=[q_return,[double(H(j)); queue(2,i)]];
                                    queue(:,i)=[];
                                end
                            end
                        end
                    end
                end
            end
        end
        patient_visitor = hospital_wait;
        hospital_wait=zeros(3,3,3);
    end
    time_stamp=[queue(2:2:end),q_return(2:2:end),queue_1(2:2:end)];
    queue=[];
    queue=[q_return queue_1];
    q_return=[];
    queue_1=[];
    hospital_wait_all=[hospital_wait_all  patient_visitor];
    patient_visitor=zeros(3,3,3);
    if n1>3
        u1=0.9;
    else
        u1=1.1;
    end
    if sum(hospital_wait_all(1,3*T+1:1:3*(T+1),1))==1
        price(1)=price(1)+2;
        time(1)=time(1)+2;
        costA=costA+2;
        timeA=timeA+2;
    elseif sum(hospital_wait_all(2,3*T+1:1:3*(T+1),1))==1
        price(1)=price(1)+5;
        time(1)=time(1)+6;
        costB=costB+5;
        timeB=timeB+6;
    elseif sum(hospital_wait_all(3,3*T+1:1:3*(T+1),1))==1
        price(1)=price(1)+7;
        time(1)=time(1)+4;
        costC=costC+7;
        timeC=timeC+4;
    end
    if  sum(hospital_wait_all(1,3*T+1:1:3*(T+1),2))==1
        kk(2)=0;
        price(2)=price(2)+7;
        time(2)=time(2)+10;
        costA=costA+7;
        timeA=timeA+10;
    elseif sum(hospital_wait_all(2,3*T+1:1:3*(T+1),2))==1
        kk(2)=2;
        price(2)=price(2)+2;
        time(2)=time(2)+5;
        costB=costB+2;
        timeB=timeB+5;
    elseif sum(hospital_wait_all(3,3*T+1:1:3*(T+1),2))==1
        kk(2)=0;
        price(2)=price(2)+5;
        time(2)=time(2)+15;
        costC=costC+5;
        timeC=timeC+15;
    else
        kk(2)=2;
    end
    if sum(hospital_wait_all(1,3*T+1:1:3*(T+1),3))==1
        if tt==0
            kk(3)=0;
            price(3)=price(3)+5;
            time(3)=time(3)+30;
            costA=costA+5;
            timeA=timeA+30;
            tt=22;
        elseif tt~=0
            tt=max(0,tt-8);
            if tt==0
                kk(3)=3;
            end
        end
    elseif sum(hospital_wait_all(2,3*T+1:1:3*(T+1),3))==1
        if ttt==0
            kk(3)=0;
            price(3)=price(3)+7;
            time(3)=time(3)+20;
            costB=costB+7;
            timeB=timeB+20;
            ttt=12;
        elseif ttt~=0
            ttt=max(0,ttt-8);
            if ttt==0
                kk(3)=3;
            end
        end
    elseif sum(hospital_wait_all(3,3*T+1:1:3*(T+1),3))==1
        kk(3)=0;
        price(3)=price(3)+2;
        time(3)=time(3)+10;
        costC=costC+2;
        timeC=timeC+10;
    else
        kk(3)=3;
    end
end
for BB=0:1:TED-1
hospital_wait_all(:,3*BB+1:1:3*(BB+1),:)*(8*(BB-1));
end
MM=1;
TA=sum(hospital_wait_all(1,:,1))+sum(hospital_wait_all(1,:,2))+sum(hospital_wait_all(1,:,3));
TB=sum(hospital_wait_all(2,:,1))+sum(hospital_wait_all(2,:,2))+sum(hospital_wait_all(2,:,3));
TC=sum(hospital_wait_all(3,:,1))+sum(hospital_wait_all(3,:,2))+sum(hospital_wait_all(3,:,3));
while(MM<=length(hospital_wait_all))
    if sum(sum(hospital_wait_all(:,MM,:)))==0
        hospital_wait_all(:,MM,:)=[];
        MM=MM-1;
    end
    MM=MM+1;
end
mean_time2=(sum(sum(sum(8*hospital_wait_all)))-sum(time_stamp))/length(time_stamp);
fprintf('\t\t\tmean_time   \tmean_cost\n')
fprintf('\t\t A \t %f \t   %f \n',timeA/TA,costA/TA)
fprintf('\t\t B \t %f \t   %f \n',timeB/TB,costB/TB)
fprintf('\t\t C \t %f \t   %f \n\n\n', timeC/TC,costC/TC)
% price=price./TED;
% time=time./TED;
% fprintf('\n\na\t%f\t%f\nb\t%f\t%f\nc\t%f\t%f\n',price(1),time(1),price(2),time(2),price(3),time(3));







