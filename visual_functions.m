classdef visual_functions
    methods
        function  outcome_loc=animate_data(obj,file_name,x,ce_save,fig_nb,xlab,ylab,title_)
            flag = 1;
            outcome_loc=flag;
            close all;
            fs = 16;

            %% Animation

            cd('Videos')
            v = VideoWriter(file_name);
            v.FrameRate = 8;
            open(v);

            temp=size(ce_save);
            Nn=temp(1);
            NT=temp(2);

            for k = 1:NT

                if ce_save(:,k)~=ones(Nn,1)

                    figure(fig_nb);
                    %set(gcf,'Position');%,[291 85 875 578]);


                    cla;
                    plot(x,ce_save(:,k),'r-','LineWidth',2);
                    %xlim([1 Nn]);
                    %ylim([min(min(ce_save)) max(max(ce_save))]);
                %     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
                %     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
                %     ylabel('$$\phi_{s,n}$$','interpreter','latex')
                    ylabel(ylab,'FontSize',fs);
                %     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                    xlabel(xlab,'FontSize',fs);
                %     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                    %set(gcf,'FontSize', fs)
                    %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                    title(title_,'fontsize',fs);
                    set(gca,'FontSize', fs)

                    tit_time = sprintf('\b Time : %3.0f',k);
                    title(tit_time,'fontsize',fs);

                    % Save Frame into video

                    frame = getframe(gcf);
                    writeVideo(v,frame);
                    pause(0.01);
                end
            end

            %% Close Video

            close(v);
            cd("..")

        end
        
        function  outcome_loc=animate_data_solid(obj,file_name,x,csn_save,csp_save,fig_nb,xlab,ylab,title_,nnb,snb,pnb)
            flag = 1;
            outcome_loc=flag;
            close all;
            fs = 16;

            %% Animation

            cd('Videos')
            v = VideoWriter(file_name);
            v.FrameRate = 8;
            open(v);

            temp=size(csn_save);
            Nn=temp(1);
            NT=temp(2);

            for k = 1:NT

                if csn_save(:,k)~=ones(Nn,1)

                    figure(fig_nb);
                    %set(gcf,'Position');%,[291 85 875 578]);

                    cla;
                    
                    plot(x(1:nnb),csn_save(:,k),'LineWidth',2);
                    hold on
                    plot(x(nnb+snb+1:nnb+snb+pnb),csp_save(:,k),'LineWidth',2);
                    %xlim([1 Nn]);
                    %ylim([min(min(csn_save)) max(max(csn_save))]);
                    %ylim([min(min(phi_s_n)) max(max(phi_s_n))])
                    %ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
                    %ylabel('$$\phi_{s,n}$$','interpreter','latex')
                    ylabel(ylab,'FontSize',fs);
                    %xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                    xlabel(xlab,'FontSize',fs);
                    %legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                    %set(gcf,'FontSize', fs)
                    %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                    title(title_,'fontsize',fs);

                    set(gca,'FontSize', fs)

                    tit_time = sprintf('\b Time : %3.0f',k);
                    title(tit_time,'fontsize',fs);

                    % Save Frame into video

                    frame = getframe(gcf);
                    writeVideo(v,frame);
                    pause(0.01);
                end
            end

            %% Close Video

            close(v);
            cd("..")

        end

        function  outcome_loc=animate_complete_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array,voltage)
            global hist
            flag = 1;
            outcome_loc=flag;
            close all;
            fs = 16;

            %% Animation

            cd('Videos')
            disp("debug vis")
            disp(file_name)
            v = VideoWriter(file_name);
            v.FrameRate = 8;
            open(v);

            temp=size(hist.ce);
            Nn=temp(1);
            NT=temp(2);
            temp=size(hist.cse);
            Nn2=temp(1);

            historic_limit=0;

            for k = 1:NT
                figure(fig_nb);
                set(gcf,'Position',[50 50 1800 1000]);
                tit_time = sprintf('\b ite : %1.0f',k);

                time_array_=cat(2,[0],time_array);
                voltage_= cat(2,[0],voltage);

                % Create a uicontrol of type "text"
                mTextBox = uicontrol('style','text');
                titleStr = sprintf('Time : %2.0f sec | Voltage : %1.2f V | iteration : %1.f',time_array_(k),voltage_(k),k);
                set(mTextBox,'String',titleStr,'FontSize',fs,'Position',[650 950 500 40]);

                %title(tit_time,'fontsize',fs);
                
                if hist.ce(:,k)~=zeros(Nn,1)

                    subplot(3,4,[1,2])

                    cla;
                    plot(x,hist.ce(:,k),'r-','LineWidth',2);
                    %xlim([1 Nn]);


                    if historic_limit==1
                        ylim([min(min(hist.ce)) max(max(hist.ce))]);
                    else
                        ylim([min(hist.ce(:,k)) max(min(hist.ce(:,k))+0.0000000001,max(hist.ce(:,k)))]);
                    end
                %     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
                %     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
                %     ylabel('$$\phi_{s,n}$$','interpreter','latex')
                    ylabel('Concentration','FontSize',fs);
                %     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                    xlabel('x','FontSize',fs);
                %     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                    %set(gcf,'FontSize', fs)
                    %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                    title('Li concentration in electrolyte','fontsize',fs);
                    set(gca,'FontSize', fs)

                    % Save Frame into video

                    frame = getframe(gcf);
                    writeVideo(v,frame);
                    
                end
                
                if hist.cse(:,k)~=zeros(Nn2,1)
                    %set(gcf,'Position');%,[291 85 875 578]);
                    subplot(3,4,[3,4])
                    
                    cla;
                    %plot(x,hist.cse(:,k),'r-','LineWidth',2);
                    plot(x(1:nnb),hist.cse(1:nnb,k),'LineWidth',2);
                    hold on
                    plot(x(nnb+snb+1:nnb+snb+pnb),hist.cse(nnb+1:nnb+pnb,k),'LineWidth',2);
                    %xlim([1 Nn]);
                    if historic_limit==1
                        ylim([min(min(hist.cse)) max(max(hist.cse))]);
                    else
                        ylim([min(hist.cse(:,k)) max(hist.cse(:,k))]);
                    end
                %     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
                %     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
                %     ylabel('$$\phi_{s,n}$$','interpreter','latex')
                    ylabel('Concentration','FontSize',fs);
                %     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                    xlabel('x','FontSize',fs);
                %     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                    %set(gcf,'FontSize', fs)
                    %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                    title('Li concentration at solid particles surface','fontsize',fs);
                    set(gca,'FontSize', fs)

                    % Save Frame into video

                    frame = getframe(gcf);
                    writeVideo(v,frame);
                    
                end
                
                if hist.pe(:,k)~=zeros(Nn,1)
                    %set(gcf,'Position');%,[291 85 875 578]);
                    subplot(3,4,[5,6])
                    
                    cla;
                    plot(x,hist.pe(:,k),'LineWidth',2);
                    %xlim([1 Nn]);
                    if historic_limit==1
                        ylim([min(min(hist.pe(:,2:NT))) max(max(hist.pe(:,2:NT)))]);
                    else
                        ylim([min(hist.pe(:,k)) max(max(hist.pe(:,k)),min(hist.pe(:,k))+0.00000000001)]);
                    end
                %     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
                %     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
                %     ylabel('$$\phi_{s,n}$$','interpreter','latex')
                    ylabel('Potential [V]','FontSize',fs);
                %     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                    xlabel('x','FontSize',fs);
                %     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                    %set(gcf,'FontSize', fs)
                    %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                    title('Potential in the electrolyte','fontsize',fs);
                    set(gca,'FontSize', fs)


                    % Save Frame into video

                    frame = getframe(gcf);
                    writeVideo(v,frame);
                    
                end
            

                if hist.ps(:,k)~=zeros(Nn2,1)
                    %set(gcf,'Position');%,[291 85 875 578]);
                    subplot(3,4,[7,8])
                    
                    cla;
                    %plot(x,hist.cse(:,k),'r-','LineWidth',2);
                    plot(x(1:nnb),hist.ps(1:nnb,k),'LineWidth',2);
                    hold on
                    plot(x(nnb+snb+1:nnb+snb+pnb),hist.ps(nnb+1:nnb+pnb,k),'LineWidth',2);
                    %xlim([1 Nn]);
                    if historic_limit==1
                        ylim([min(min(hist.ps)) max(max(hist.ps))]);
                    else
                        ylim([min(hist.ps(:,k)) max(hist.ps(:,k))]);
                    end
                %     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
                %     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
                %     ylabel('$$\phi_{s,n}$$','interpreter','latex')
                    ylabel('Potential [V]','FontSize',fs);
                %     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                    xlabel('x','FontSize',fs);
                %     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                    %set(gcf,'FontSize', fs)
                    %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                    title('Potential in the solid','fontsize',fs);
                    set(gca,'FontSize', fs)

                    % Save Frame into video

                    frame = getframe(gcf);
                    writeVideo(v,frame);
                    
                end
                
                
                %set(gcf,'Position');%,[291 85 875 578]);
                subplot(3,4,[9,10])

                cla;
                plot(x,hist.j(:,k),'LineWidth',2);
                %xlim([1 Nn]);
                if historic_limit==1
                    ylim([min(min(hist.j(:,2:NT))) max(max(hist.j(:,2:NT)))]);
                else
                    ylim([min(hist.j(:,k)) max(max(hist.j(:,k)),min(hist.j(:,k))+0.000000001)]);
                end
             
            %     ylim([min(min(phi_s_n)) max(max(phi_s_n))])
            %     ylabel('$$c_{s,n}(x,t)$$','interpreter','latex','FontSize',fs)
            %     ylabel('$$\phi_{s,n}$$','interpreter','latex')
                ylabel('j','FontSize',fs);
            %     xlabel('Space across cell, x [node no.]','interpreter','latex','FontSize',fs)
                xlabel('x','FontSize',fs);
            %     legend({'$$c_{ss,n}(x,t)$$';'$$c_{avg,n}(x,t)$$'},'interpreter','latex','FontSize',fs,0);
                %set(gcf,'FontSize', fs)
                %set(gcf,'Position',[0.1300 0.7093 0.31 0.2157]);
                title('rate of positive charge flowing','fontsize',fs);
                set(gca,'FontSize', fs)

                % Save Frame into video

                frame = getframe(gcf);
                writeVideo(v,frame);
                    
                
                
                
            end
            
            %% Close Video

            close(v);
            cd('..')

        end
        
        function outcome_loc= plot_data(obj,x,y,xlab,ylab,title_,fig_nb,file_name)
            cd('Graphs')
            figure(fig_nb);
            clf;
            fs = 16;
            %subplot(4,1,1);
            plot(x,y,'LineWidth',2);
            ylabel(ylab,'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$$I(t)$$'},'interpreter','latex');
            title(title_,'fontsize',fs);
            xlabel(xlab,'FontSize',fs);
            saveas(gcf,file_name);
            cd("..")
        end


        function outcome_loc= plot_data_solid(obj,x,y,xlab,ylab,title_,fig_nb,file_name,nnb,snb,pnb)
            cd('Graphs')
            figure(fig_nb);
            clf;
            fs = 16;
            %subplot(4,1,1);
            plot(x(1:nnb),y(1:nnb),'LineWidth',2);
            hold on
            plot(x(nnb+snb+1:nnb+snb+pnb),y(nnb+1:nnb+pnb),'LineWidth',2);
            ylabel(ylab,'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$$I(t)$$'},'interpreter','latex');
            title(title_,'fontsize',fs);
            xlabel(xlab,'FontSize',fs);
            saveas(gcf,file_name);
            cd("..")
        end
        
        function  outcome_loc=plot_complete_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array)
            global fv
            global hist
            flag = 1;
            outcome_loc=flag;
            close all;
            fs = 16;

            cd('Graphs')
            
            temp=size(fv.ce);
            Nn=temp(1);
            temp=size(fv.cse);
            Nn2=temp(1);

            
            figure(fig_nb);
            set(gcf,'Position',[50 50 1800 1000]);     
            
            subplot(5,5,[1,2])

            cla;
            plot(x,fv.ce,'r-','LineWidth',2);

            if min(fv.ce) == max(fv.ce)
                ylim([min(fv.ce) min(fv.ce)+min(fv.ce)*0.00001]);
            else
                ylim([min(fv.ce) max(fv.ce)]);
            end

            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in electrolyte','fontsize',fs);
            set(gca,'FontSize', fs)
            
        

            subplot(5,5,[4,5])

            cla;
            plot(x,fv.pe,'LineWidth',2);
            ylim([min(fv.pe(:)) max(max(fv.pe(:)),min(fv.pe(:))+abs(min(fv.pe(:)))*0.0001)]);
        
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            
            title('Potential in the electrolyte','fontsize',fs);
            set(gca,'FontSize', fs)



            subplot(5,5,[6,7])

            cla;
            plot(x(1:nnb),fv.cse(1:nnb),'LineWidth',2);
            ylim([min(fv.cse(1:nnb)) max(max(fv.cse(1:nnb)),min(fv.cse(1:nnb))+0.0001)]);
        
            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Li concentration at solid particles surface','fontsize',fs);
            set(gca,'FontSize', fs)



            subplot(5,5,[9,10])

            cla;
            plot(x(nnb+snb+1:nnb+snb+pnb),fv.cse(nnb+1:nnb+pnb),'LineWidth',2);
            ylim([min(fv.cse(nnb+1:nnb+pnb)) max(max(fv.cse(nnb+1:nnb+pnb)),min(fv.cse(nnb+1:nnb+pnb))+0.0001)]);
        
            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Li concentration at solid particles surface','fontsize',fs);
            set(gca,'FontSize', fs)
        

            
            subplot(5,5,[11,12])

            cla;
            plot(x(1:nnb),fv.ps(1:nnb),'LineWidth',2);
            
            if max(fv.ps)==0
                ylim([min(fv.ps(1:nnb)) 0.00000000001]);
            else
                ylim([min(fv.ps(1:nnb)) max(max(fv.ps(1:nnb)),min(fv.ps(1:nnb))+0.00000000001)]);
            end
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Potential in the solid','fontsize',fs);
            set(gca,'FontSize', fs)        

            
            subplot(5,5,[14,15])

            cla;
            plot(x(1:nnb),fv.ps(1:nnb),'LineWidth',2);
            hold on
            plot(x(nnb+snb+1:nnb+snb+pnb),fv.ps(nnb+1:nnb+pnb),'LineWidth',2);
            if max(fv.ps)==0
                ylim([min(fv.ps(nnb+1:nnb+pnb)) 0.00000000001]);
            else
                %disp([min(fv.ps(nnb+1:nnb+pnb)) max(max(fv.ps(nnb+1:nnb+pnb)),min(fv.ps(nnb+1:nnb+pnb))+0.00000000001)])
                ylim([min(fv.ps(nnb+1:nnb+pnb)) max(max(fv.ps(nnb+1:nnb+pnb)),min(fv.ps(nnb+1:nnb+pnb))+0.00000000001)]);
            end
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Potential in the solid','fontsize',fs);
            set(gca,'FontSize', fs)




            subplot(5,5,[16,17])

            cla;
            plot(x(1:nnb),fv.j(1:nnb),'LineWidth',2);
            if max(fv.j(1:nnb))==0
                ylim([min(fv.j(1:nnb)) 0.00000000001]);
            else
                ylim([min(fv.j(1:nnb)) max(max(fv.j(1:nnb)),min(fv.j(1:nnb))+abs(min(fv.j(1:nnb)))*0.00001)]);
            end
            ylabel('j','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('rate of positive charge flowing','fontsize',fs);
            set(gca,'FontSize', fs)


            subplot(5,5,[19,20])

            cla;
            plot(x(nnb+snb+1:nnb+snb+pnb),fv.j(nnb+snb+1:nnb+snb+pnb),'LineWidth',2);
            if max(fv.j(nnb+snb+1:nnb+snb+pnb))==0
                ylim([min(fv.j(nnb+snb+1:nnb+snb+pnb)) 0.00000000001]);
            else
                ylim([min(fv.j(nnb+snb+1:nnb+snb+pnb)) max(max(fv.j(nnb+snb+1:nnb+snb+pnb)),min(fv.j(nnb+snb+1:nnb+snb+pnb))+abs(min(fv.j(nnb+snb+1:nnb+snb+pnb)))*0.0001)]);
            end
            ylabel('j','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('rate of positive charge flowing','fontsize',fs);
            set(gca,'FontSize', fs)


            subplot(5,5,[21,22])

            cla;
            plot(time_array,hist.V,'LineWidth',2);
            if max(hist.V(:))==0
                ylim([min(hist.V(:)) 0.0000001]);
            else
                ylim([min(hist.V(:)) max(max(hist.V(:)),min(hist.V(:))+abs(min(hist.V(:)))*0.0001)]);
            end
            ylabel('Voltage [V]','FontSize',fs);
            xlabel('time [sec]','FontSize',fs);
        
            title('Cell voltage over time','fontsize',fs);
            set(gca,'FontSize', fs)



            saveas(gcf,file_name);


            cd('..')

        end


        function outcome_loc= plot_resuduals_newt(obj,file_name)
            global hist
            global sol
            cd('Graphs')
            figure(234234234);
            set(gcf,'Position',[50 50 1800 1000]);     
            clf;
            fs = 16;

            subplot(3,4,[1,2]);

            
            cla;
            x=1:hist.newt_it_number(1,sol.nb_steps);
            y=hist.residuals(1,1:hist.newt_it_number(1,sol.nb_steps));
            graphic = semilogy(x,y,'LineWidth',2);
            ylabel("residual coupled",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[3,4]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.nb_steps);
            y=hist.residuals(2,1:hist.newt_it_number(1,sol.nb_steps));
            graphic = semilogy(x,y,'LineWidth',2);
            ylabel("residual ps",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);

            subplot(3,4,[5,6]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.nb_steps);
            y=hist.residuals(3,1:hist.newt_it_number(1,sol.nb_steps));
            graphic = semilogy(x,y,'LineWidth',2);
            ylabel("residual pe",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[7,8]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.nb_steps);
            y=hist.residuals(4,1:hist.newt_it_number(1,sol.nb_steps));
            graphic = semilogy(x,y,'LineWidth',2);
            ylabel("residual ce",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[9,10]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.nb_steps);
            y=hist.residuals(5,1:hist.newt_it_number(1,sol.nb_steps));
            graphic = semilogy(x,y,'LineWidth',2);
            ylabel("residual csn",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[11,12]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.nb_steps);
            y=hist.residuals(6,1:hist.newt_it_number(1,sol.nb_steps));
            graphic = semilogy(x,y,'LineWidth',2);
            ylabel("residual csp",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);

            saveas(gcf,file_name);
            cd("..")
        end

        function outcome_loc= plot_resuduals_DFN(obj,file_name)
            global hist
            global sol
            cd('Graphs')
            figure(234234235);
            set(gcf,'Position',[50 50 1800 1000]);     
            clf;
            fs = 16;

            subplot(4,4,[1,2]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_time(1,:);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual coupled",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(4,4,[3,4]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_time(2,:);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual ps",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(4,4,[5,6]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_time(3,:);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual pe",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(4,4,[7,8]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_time(4,:);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual ce",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(4,4,[9,10]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_time(5,:);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual csn",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(4,4,[11,12]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_time(6,:);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual csp",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);


            saveas(gcf,file_name);
            cd("..")
        end

        function outcome_loc= plot_resuduals_diff(obj,file_name)
            global hist
            global sol
            cd('Graphs')
            figure(234234236);
            set(gcf,'Position',[50 50 1800 1000]);     
            clf;
            fs = 16;

            subplot(4,4,[1,2]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_diff(1,:);
            plot(x,y,'LineWidth',2);
            ylabel("residual coupled",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[3,4]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_diff(2,:);
            plot(x,y,'LineWidth',2);
            ylabel("residual ps",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[5,6]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_diff(3,:);
            plot(x,y,'LineWidth',2);
            ylabel("residual pe",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[7,8]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_diff(4,:);
            plot(x,y,'LineWidth',2);
            ylabel("residual ce",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[9,10]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_diff(5,:);
            plot(x,y,'LineWidth',2);
            ylabel("residual csn",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[11,12]);
            cla;
            x=1:sol.nb_steps;
            y=hist.residuals_diff(6,:);
            plot(x,y,'LineWidth',2);
            ylabel("residual csp",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);


            subplot(4,4,[13,14]);
            cla;
            x=1:sol.nb_steps;
            y=hist.newt_it_number(1,:);
            plot(x,y,'LineWidth',2);
            ylabel("Newton solver iterations",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("iteration",'FontSize',fs);


            saveas(gcf,file_name);
            cd("..")
        end


        function  outcome_loc=plot_solid_concentration_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array)
            global fv
            global hist
            global sol
            flag = 1;
            outcome_loc=flag;
            fs = 16;

            cd('Graphs')
            
            temp=size(hist.csn);
            Nn=temp(2);
            temp=size(hist.csp);
            Nn2=temp(2);

            
            figure(234234237);
            set(gcf,'Position',[50 50 1800 1000]);     
            
            subplot(2,2,1)
            x=1:sol.nb_cell_n*(sol.part_nb_cell+1);
            y=reshape(fv.csn,sol.nb_cell_n*(sol.part_nb_cell+1),1);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+min(y)*0.00001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csn','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in neg. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
            

            subplot(2,2,2)
            x=1:sol.nb_cell_p*(sol.part_nb_cell+1);
            y=reshape(fv.csp,sol.nb_cell_p*(sol.part_nb_cell+1),1);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+min(y)*0.00001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csp','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        

            subplot(2,2,3)
            x=1:sol.nb_cell_n*(sol.part_nb_cell+1);
            y=hist.csn(:,Nn-1);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csn -','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        

            subplot(2,2,4)
            x=1:sol.nb_cell_p*(sol.part_nb_cell+1);
            y=hist.csp(:,Nn2-1);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csp -','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        
            


            saveas(gcf,file_name);


            cd('..')

        end


        function  outcome_loc=plot_solid_concentration_singlePart_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array)
            global fv
            global hist
            global sol
            flag = 1;
            outcome_loc=flag;
            fs = 16;

            cd('Graphs')

            temp=size(hist.csn);
            Nn=temp(2);
            temp=size(hist.csp);
            Nn2=temp(2);

            
            figure(234234247);
            set(gcf,'Position',[50 50 1800 1000]);     
            
            subplot(2,2,1)
            x=sol.part_coord_n;
            y=hist.csn(1:sol.part_nb_cell+1,Nn);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csn1','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in neg. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
            

            subplot(2,2,2)
            x=sol.part_coord_n;
            y=hist.csn(sol.part_nb_cell+2:2*sol.part_nb_cell+2,Nn);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csn2','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        

            subplot(2,2,3)
            x=sol.part_coord_p;
            y=hist.csp(1:sol.part_nb_cell+1,Nn2);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csp1','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        

            subplot(2,2,4)
            x=sol.part_coord_p;
            y=hist.csp(sol.part_nb_cell+2:2*sol.part_nb_cell+2,Nn2);

            cla;
            plot(x,y,'r-','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csp2','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        
            


            saveas(gcf,file_name);


            cd('..')
        end



    end
end
