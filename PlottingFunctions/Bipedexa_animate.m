function Bipedexa_animate(output,savename,varargin)
% Makes a 4-second animation of the result
% Requires patchline:
% http://www.mathworks.com/matlabcentral/fileexchange/36953-patchline/content/patchline.m


p = inputParser;
addRequired(p,'output',@isstruct)
addRequired(p,'savename',@isstr)

addParameter(p,'alpha',true,@islogical)
addParameter(p,'showtext',[0 0 0],@isvector)
addParameter(p,'threshold',0.05,@(x) isscalar(x) && isnumeric(x))
addParameter(p,'BodyType','DistPM',@isstr)
addParameter(p,'txt',true)

parse(p,output,savename,varargin{:})

alpha = p.Results.alpha;
showtext = p.Results.showtext;
threshold = p.Results.threshold;
bodytype = p.Results.BodyType;
txt = p.Results.txt;

filename = [savename,'_animate'];

auxdata = output.result.setup.auxdata;
T = auxdata.T;
tq = linspace(0,T,120)'; % query time points; 120 frames


D = auxdata.D;
d = auxdata.d;

t = output.result.interpsolution.phase.time;
X = output.result.interpsolution.phase.state;

x = X(:,1);
y = X(:,2);
theta = X(:,5);
F = X(:,7:9);

Fq = interp1(t,F,tq,'pchip');
xq = interp1(t,x,tq,'pchip');
yq = interp1(t,y,tq,'pchip');
thetaq = interp1(t,theta,tq,'pchip');

% Make new inertial reference frame at average horizontal speed
x0 = tq*D/T;

% Get new footfall locations
f_tr = -x0;
f_ld = D-x0;
f_rf = d-x0;

Fmax = max(max(F));
lmax = auxdata.lmax;
r = auxdata.r;
lt = lmax; % for now, make the torso length equal to leg length. This allows easy plotting, but the perceived MOI will be wrong.
m = [r/lt,1 - r/lt];
F = abs(Fq./Fmax);
threshold = threshold/Fmax; % as force is normalized to Fmax, normalize the threshold too
Fr = auxdata.Fr;
Uh = D*sqrt(Fr/lmax);
c = auxdata.c;
I = auxdata.I;



rs = (lt-r)*[cos(thetaq), sin(thetaq)];
rh = -r*[cos(thetaq), sin(thetaq)];

xsh = [rs(:,1)+xq-x0,rh(:,1)+xq-x0]; % x position of shoulders and hips
ysh = [rs(:,2)+yq,rh(:,2)+yq]; % y position of shoulders and hips


% Set axis limits
yl = [0, 1.5*(lmax+lt)];
xl = [-2.5 2.5];%[-1.5,1.5];



close all;
figure('Position', [440 378 560 2*lmax*560/3],'color','w')
axes('Position',[.1 .1 .8 .8])
hold on;

switch upper(bodytype)
    case {'EVENPM','EVENPOINTMASS'}
    rgyr = sqrt(I)*[cos(thetaq), sin(thetaq)]/2; % radius of gyration, relative to half body length.
    xgyr = [-rgyr(:,1)+xq-x0,+rgyr(:,1)+xq-x0];
    ygyr = [-rgyr(:,2)+yq,rgyr(:,2)+yq];
    mrksz = [10 10]; % marker sizes are equal
    mrkst = 's';
    case {'PMLIMBS','POINTMASSLIMBS','DISTPM','DISTRIBUTEDPOINTMASS'}
        xgyr = xsh;
        ygyr = ysh;
        mrksz = 25*m;
        mrkst = 'o';
    case {'PM','POINTMASS'}
        xgyr = (xq-x0)*[1 1];
        ygyr = [yq,yq];
        mrksz = 18*[1 1];
        mrkst = 'o';
end

Work = output.result.solution.phase.integral(1);
slackPen = output.result.solution.phase.integral(2);



writerObj = VideoWriter(filename,'MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

[txt1,txt2,txt3] = deal('');
if ischar(showtext)
    switch upper(showtext)
        case 'OUTSIDE'
            if isnan(txt)
                txt1 = ['$U'' = ', num2str(Uh,'%.2f'),',\quad','\hat{I} = ',num2str(Ihat,2),'$'];
            else
                txt1 = txt;
            end
            xy_txt = [0.5,1.05;
                      0.5,0.9;
                      0.9,0.9];
            alignment_txt = {'center','bottom'};
            fontsize_txt = 14;
        case 'OUTSIDE+INSIDE'
            if iscell(txt)
                txt1 = txt{1};
                txt2 = txt{2};
            elseif isnan(txt)
                txt1 = ['$U'' = ', num2str(Uh,'%.2f'),',\quad','\hat{I} = ',num2str(Ihat,2),'$'];
                txt2 = ['$D'' = $',num2str(D),char(10)...
                        '$l''_{Fmax} = $', num2str(lmaxF),char(10)];
            end
            xy_txt = [0.5,1.05;
                      0.5,0.9;
                      0.9,0.9];
            alignment_txt = {'center','bottom'};
            fontsize_txt = 14;
    end
else
    showtext = logical(showtext);
    xy_txt = [0.1,0.9;
              0.25,0.9;
              0.9,0.9];
    alignment_txt = {'left','top'};
    fontsize_txt = 12;
    if showtext(1)
        % This info to be plotted at top-left
        txt1 = ['$U''_H$ = ', num2str(Uh,2),newline,...
            '$D'' = $',num2str(D),newline...
            '$l''_{Fmax} = $', num2str(lmaxF),newline];
    end
    if showtext(2)
        % This info to be plotted to the right of the above
        c_text = ['[',regexprep(num2str(c), ' *', ' '),']'];
        txt2 = ['nlpinfo: ',num2str(output.result.nlpinfo),newline,...
            'c = ', c_text,...
            newline,' CoT: Work: ', num2str(Work/D,4),' $\dot{F}^2$: ', num2str(Fdot/D,4), ' slackPen: ' num2str(slackPen,4)];
    end
    if showtext(3)
        txt3 = ['Threshold = $F_{peak}\times$', num2str(threshold)];
    end
end
lw = 2;
for i = 1:length(tq)
    cla;

    % Parameters
    text(xy_txt(1,1),xy_txt(1,2),txt1,'interpreter','latex','units','normalized','horizontalalignment',alignment_txt{1},'verticalalignment',alignment_txt{2},'fontsize',fontsize_txt);
    text(xy_txt(2,1),xy_txt(2,2),txt2,'interpreter','latex','units','normalized','horizontalalignment',alignment_txt{1},'verticalalignment','top','fontsize',fontsize_txt);
    text(xy_txt(3,1),xy_txt(3,2),txt3,...
        'units','normalized','verticalalignment','top',...
        'horizontalalignment','right','interpreter','latex','fontsize',fontsize_txt);
    % limbs
    if alpha
        patchline([f_tr(i),xsh(i,2)],[0,ysh(i,2)],'linestyle','--','edgecolor','b','edgealpha',F(i,1),'linewidth',lw);
        patchline([f_ld(i),xsh(i,2)],[0,ysh(i,2)],'linestyle','--','edgecolor','b','edgealpha',F(i,2),'linewidth',lw);
        patchline([f_rf(i),xsh(i,2)],[0,ysh(i,2)],'linestyle','-','edgecolor','r','edgealpha',F(i,3),'linewidth',lw);
    else
        if Fq(i,1) > threshold
            patchline([f_tr(i),xsh(i,2)],[0,ysh(i,2)],'linestyle','--','edgecolor','b','linewidth',lw)
        end
        if Fq(i,2) > threshold
            patchline([LF1(i),xsh(i,1)],[0,ysh(i,1)],'linestyle','--','edgecolor','b','linewidth',lw)
        end
        if Fq(i,3) > threshold
            patchline([LF2(i),xsh(i,1)],[0,ysh(i,1)],'linestyle','--','edgecolor','b','linewidth',lw)
        end
    end
    
    % Torso
    plot(xsh(i,:),ysh(i,:),'k:','linewidth',lw-0.5)
    plot(xgyr(i,:),ygyr(i,:),'k-','linewidth',lw)
    plot(xgyr(i,1),ygyr(i,1),'ks','linewidth',lw,'markersize',mrksz(1),'marker',mrkst,'markerfacecolor','w')
    plot(xgyr(i,2),ygyr(i,2),'ks','linewidth',lw,'markersize',mrksz(2),'marker',mrkst,'markerfacecolor','w')
    plot(xq(i)-x0(i),yq(i),'kx','markersize',10)
    
 
    % ground markers
    plot(f_tr(i),0,'bo')
    
    plot(f_ld(i),0,'bs')
    plot(f_rf(i),0,'ro')
    
    plot(xl,[0,0],'k:')
    axis equal
    ylim(yl);
    xlim(xl);
    box on
    xlabel('x/l_b')
    ylabel('y/l_b')
    set(gca,'linewidth',1)
    
    drawnow
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end
close(writerObj)
end