

classdef visual_functions
    methods
        function  outcome_loc=animate_data(obj,file_name,x,ce_save,fig_nb,xlab,ylab,title_)
            global deb
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
                    
                    grid on
                    grid minor
                    
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
            global deb
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
                    %set(gcf,'Position'); %,[291 85 875 578]);

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

                    grid on
                    grid minor

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
            global deb
            flag = 1;
            outcome_loc=flag;
            close all;
            fs = 16;

            %% Animation

            cd('Videos')
            
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
        
        function outcome_loc= plot_data(obj,x,y,xlab,ylab,title_,fig_nb,file_name,compare_file_name,expe_file_name)
            global deb
            plot_comp=0;
            plot_expe=0;
            if not(compare_file_name=="")
                plot_comp=1;
                Tcomp = table2array(readtable(compare_file_name));
            end
            if not(expe_file_name=="")
                plot_expe=1;
                Texpe = table2array(readtable(expe_file_name));
            end

            figure('visible','off');
            clf;
            fs = 16;
            %subplot(4,1,1);
            if plot_comp==1
                plot(Tcomp(:,1),Tcomp(:,2),'LineWidth',2);
                hold on
            end
            if plot_expe==1
                plot(Texpe(:,1),Texpe(:,2),'LineWidth',2);
                hold on
            end
            plot(x,y,'LineWidth',2);

            if plot_comp==1
                legend('Other model','expe','DFN DianNiu')
            else
                legend('DFN DianNiu')                
            end

            ylabel(ylab,'FontSize', fs);
            set(gca,'FontSize', fs);
            title(title_,'fontsize',fs);
            xlabel(xlab,'FontSize',fs);
            grid on
            grid minor
            saveas(gcf,deb.graph_folder_name+file_name);
        end


        function outcome_loc= plot_data_solid(obj,x,y,xlab,ylab,title_,fig_nb,file_name,nnb,snb,pnb)
            global deb
            figure('visible','off');
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
            grid on
            grid minor
            saveas(gcf,deb.graph_folder_name+file_name);
        end
        
        function  outcome_loc=plot_complete_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array)
            global fv
            global hist
            global deb
            global sol
            flag = 1;
            outcome_loc=flag;
            fs = 16;

            
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
                ylim([min(fv.ce) min(fv.ce)+0.00001]);
            else
                ylim([min(fv.ce) max(fv.ce)]);
            end

            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in electrolyte','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)
            
        

            subplot(5,5,[4,5])

            cla;
            plot(x,fv.pe,'LineWidth',2);
            ylim([min(fv.pe(:)) max(max(fv.pe(:)),min(fv.pe(:))+abs(min(fv.pe(:)))*0.0001)]);
        
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            
            title('Potential in the electrolyte','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)



            subplot(5,5,[6,7])

            cla;
            plot(x(1:nnb),fv.cse(1:nnb),'LineWidth',2);
            ylim([min(fv.cse(1:nnb)) max(max(fv.cse(1:nnb)),min(fv.cse(1:nnb))+0.0001)]);
        
            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Li concentration at solid particles surface','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)



            subplot(5,5,[9,10])

            cla;
            plot(x(nnb+snb+1:nnb+snb+pnb),fv.cse(nnb+1:nnb+pnb),'LineWidth',2);
            ylim([min(fv.cse(nnb+1:nnb+pnb)) max(max(fv.cse(nnb+1:nnb+pnb)),min(fv.cse(nnb+1:nnb+pnb))+0.0001)]);
        
            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Li concentration at solid particles surface','fontsize',fs);
            grid on
            grid minor
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
            grid on
            grid minor
            set(gca,'FontSize', fs)        

            
            subplot(5,5,[14,15])

            cla;
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
            grid on
            grid minor
            set(gca,'FontSize', fs)


            
            if sum(isnan(fv.j))==0 && sum(isinf(fv.j))==0
                subplot(5,5,[16,17])

                cla;
                plot(x(1:nnb),fv.j(1:nnb),'LineWidth',2);
                if max(fv.j(1:nnb))==0
                    ylim([min(fv.j(1:nnb)) 0.00000000001]);
                else
                    ylim([min(fv.j(1:nnb)) max(max(fv.j(1:nnb)),min(fv.j(1:nnb))+abs(min(fv.j(1:nnb)))*0.00001+0.00000001)]);
                end
                ylabel('j','FontSize',fs);
                xlabel('x','FontSize',fs);
            
                title('rate of positive charge flowing','fontsize',fs);
                grid on
                grid minor
                set(gca,'FontSize', fs)

                subplot(5,5,[19,20])

                cla;
                plot(x(nnb+snb+1:nnb+snb+pnb),fv.j(nnb+snb+1:nnb+snb+pnb),'LineWidth',2);
                if max(fv.j(nnb+snb+1:nnb+snb+pnb))==0
                    ylim([min(fv.j(nnb+snb+1:nnb+snb+pnb)) 0.00000000001]);
                else
                    ylim([min(fv.j(nnb+snb+1:nnb+snb+pnb)) max(max(fv.j(nnb+snb+1:nnb+snb+pnb)),min(fv.j(nnb+snb+1:nnb+snb+pnb))+abs(min(fv.j(nnb+snb+1:nnb+snb+pnb)))*0.0001+0.00000001)]);
                end
                ylabel('j','FontSize',fs);
                xlabel('x','FontSize',fs);
            
                title('rate of positive charge flowing','fontsize',fs);
                grid on
                grid minor
                set(gca,'FontSize', fs)
            end


            subplot(5,5,[21,22])

            cla;
            plot(time_array(1:sol.time_ite),hist.V(1:sol.time_ite),'LineWidth',2);
            if max(hist.V(1:sol.time_ite))==0
                ylim([min(hist.V(1:sol.time_ite)) 0.0000001]);
            else
                ylim([min(hist.V(1:sol.time_ite)) max(max(hist.V(1:sol.time_ite)),min(hist.V(1:sol.time_ite))+abs(min(hist.V(1:sol.time_ite)))*0.0001)]);
            end
            xlim([time_array(1) time_array(sol.time_ite)]);
            ylabel('Voltage [V]','FontSize',fs);
            xlabel('time [sec]','FontSize',fs);
        
            title('Cell voltage over time','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)
            

            subplot(5,5,[24,25])

            cla;
            plot(time_array(1:sol.time_ite),hist.SOC_neg(1:sol.time_ite),'LineWidth',2);
            hold on
            plot(time_array(1:sol.time_ite),hist.SOC_pos(1:sol.time_ite),'LineWidth',2);
            ylim([0 1.0]);
            xlim([time_array(1) time_array(sol.time_ite)]);
            
            ylabel('SOC','FontSize',fs);
            xlabel('time [sec]','FontSize',fs);
        
            title('Cell State of Charge over time','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)



            saveas(gcf,deb.graph_folder_name+file_name);



        end


        function  outcome_loc=plot_complete_data_spec_ite(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array,iteration)
            global fv
            global hist
            global deb
            global sol
            flag = 1;
            outcome_loc=flag;
            fs = 16;

            
            temp=size(fv.ce);
            Nn=temp(1);
            temp=size(fv.cse);
            Nn2=temp(1);

            
            figure(fig_nb);
            set(gcf,'Position',[50 50 1800 1000]);     
            
            subplot(5,5,[1,2])

            cla;
            plot(x,hist.ce(:,iteration+1),'r-','LineWidth',2);

            if min(hist.ce(:,iteration+1)) == max(hist.ce(:,iteration+1))
                ylim([min(hist.ce(:,iteration+1)) min(hist.ce(:,iteration+1))+0.00001]);
            else
                ylim([min(hist.ce(:,iteration+1)) max(hist.ce(:,iteration+1))]);
            end

            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title("Li concentration in electrolyte ite"+num2str(iteration),'fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)
            
        

            subplot(5,5,[4,5])
            
            cla;
            plot(x,hist.pe(:,iteration+1),'LineWidth',2);
            ylim([min(hist.pe(:,iteration+1)) max(max(hist.pe(:,iteration+1)),min(hist.pe(:,iteration+1))+abs(min(hist.pe(:,iteration+1)))*0.0001)]);
        
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            
            title('Potential in the electrolyte','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)



            subplot(5,5,[6,7])

            cla;
            plot(x(1:nnb),hist.cse(1:nnb,iteration+1),'LineWidth',2);
            ylim([min(hist.cse(1:nnb,iteration+1)) max(max(hist.cse(1:nnb,iteration+1)),min(hist.cse(1:nnb,iteration+1))+0.0001)]);
        
            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Li concentration at solid particles surface','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)



            subplot(5,5,[9,10])

            cla;
            plot(x(nnb+snb+1:nnb+snb+pnb),hist.cse(nnb+1:nnb+pnb,iteration+1),'LineWidth',2);
            ylim([min(hist.cse(nnb+1:nnb+pnb,iteration+1)) max(max(hist.cse(nnb+1:nnb+pnb,iteration+1)),min(hist.cse(nnb+1:nnb+pnb,iteration+1))+0.0001)]);
        
            ylabel('Concentration','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Li concentration at solid particles surface','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)
        

            
            subplot(5,5,[11,12])

            cla;
            plot(x(1:nnb),hist.ps(1:nnb,iteration+1),'LineWidth',2);
            
            if max(hist.ps(:,iteration+1))==0
                ylim([min(hist.ps(1:nnb,iteration+1)) 0.00000000001]);
            else
                ylim([min(hist.ps(1:nnb,iteration+1)) max(max(hist.ps(1:nnb,iteration+1)),min(hist.ps(1:nnb,iteration+1))+0.00000000001)]);
            end
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Potential in the solid','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)        

            
            subplot(5,5,[14,15])

            cla;
            plot(x(nnb+snb+1:nnb+snb+pnb),hist.ps(nnb+1:nnb+pnb,iteration+1),'LineWidth',2);
            if max(hist.ps(:,iteration+1))==0
                ylim([min(hist.ps(nnb+1:nnb+pnb,iteration+1)) 0.00000000001]);
            else
                %disp([min(hist.ps(nnb+1:nnb+pnb,iteration+1)) max(max(hist.ps(nnb+1:nnb+pnb,iteration+1)),min(hist.ps(nnb+1:nnb+pnb,iteration+1))+0.00000000001)])
                ylim([min(hist.ps(nnb+1:nnb+pnb,iteration+1)) max(max(hist.ps(nnb+1:nnb+pnb,iteration+1)),min(hist.ps(nnb+1:nnb+pnb,iteration+1))+0.00000000001)]);
            end
            ylabel('Potential [V]','FontSize',fs);
            xlabel('x','FontSize',fs);
            title('Potential in the solid','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)

            sizzz=size(hist.j);

            if sizzz(2)>=iteration
                if sum(isnan(hist.j(:,iteration)))==0 && sum(isinf(hist.j(:,iteration)))==0
                    subplot(5,5,[16,17])

                    cla;
                    plot(x(1:nnb),hist.j(1:nnb,iteration),'LineWidth',2);
                    if max(hist.j(1:nnb,iteration))==0
                        ylim([min(hist.j(1:nnb,iteration)) 0.00000000001]);
                    else
                        ylim([min(hist.j(1:nnb,iteration)) max(max(hist.j(1:nnb,iteration)),min(hist.j(1:nnb,iteration))+abs(min(hist.j(1:nnb,iteration)))*0.00001+0.00000001)]);
                    end
                    ylabel('j','FontSize',fs);
                    xlabel('x','FontSize',fs);
                
                    title('rate of positive charge flowing','fontsize',fs);
                    grid on
                    grid minor
                    set(gca,'FontSize', fs)

                    subplot(5,5,[19,20])

                    cla;
                    plot(x(nnb+snb+1:nnb+snb+pnb),hist.j(nnb+snb+1:nnb+snb+pnb,iteration),'LineWidth',2);
                    if max(hist.j(nnb+snb+1:nnb+snb+pnb,iteration))==0
                        ylim([min(hist.j(nnb+snb+1:nnb+snb+pnb,iteration)) 0.00000000001]);
                    else
                        ylim([min(hist.j(nnb+snb+1:nnb+snb+pnb,iteration)) max(max(hist.j(nnb+snb+1:nnb+snb+pnb,iteration)),min(hist.j(nnb+snb+1:nnb+snb+pnb,iteration))+abs(min(hist.j(nnb+snb+1:nnb+snb+pnb,iteration)))*0.0001+0.00000001)]);
                    end
                    ylabel('j','FontSize',fs);
                    xlabel('x','FontSize',fs);
                
                    title('rate of positive charge flowing','fontsize',fs);
                    grid on
                    grid minor
                    set(gca,'FontSize', fs)
                end
            end


            subplot(5,5,[21,22])

            cla;
            plot(time_array(1:sol.time_ite),hist.V(1:sol.time_ite),'LineWidth',2);
            if max(hist.V(1:sol.time_ite))==0
                ylim([min(hist.V(1:sol.time_ite)) 0.0000001]);
            else
                ylim([min(hist.V(1:sol.time_ite)) max(max(hist.V(1:sol.time_ite)),min(hist.V(1:sol.time_ite))+abs(min(hist.V(1:sol.time_ite)))*0.0001)]);
            end
            xlim([time_array(1) time_array(sol.time_ite)]);
            ylabel('Voltage [V]','FontSize',fs);
            xlabel('time [sec]','FontSize',fs);
        
            title('Cell voltage over time','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)
            

            subplot(5,5,[24,25])

            cla;
            plot(time_array(1:sol.time_ite),hist.SOC_neg(1:sol.time_ite),'LineWidth',2);
            hold on
            plot(time_array(1:sol.time_ite),hist.SOC_pos(1:sol.time_ite),'LineWidth',2);
            ylim([min(0,min(min(hist.SOC_neg),min(hist.SOC_pos))) 1.0]);
            xlim([time_array(1) time_array(sol.time_ite)]);
            
            ylabel('SOC','FontSize',fs);
            xlabel('time [sec]','FontSize',fs);
        
            title('Cell State of Charge over time','fontsize',fs);
            grid on
            grid minor
            set(gca,'FontSize', fs)



            saveas(gcf,deb.graph_folder_name+file_name);



        end






        function ticks_array= logarithmic_ticks_generator(obj,minval,maxval)
            maxlog10=log10(maxval);
            minlog10=log10(minval);

            mintick=ceil(minlog10);
            maxtick=floor(maxlog10);
            

            temp_array=[mintick:1:maxtick];
            ticks_array=zeros(1,length(temp_array));

            for i=1:length(temp_array)
                ticks_array(i)=10^(temp_array(i));
            end
        end

        function outcome_loc= plot_resuduals_newt(obj,file_name,last_ite)
            global hist
            global sol
            global deb
            figure(234234234);
            set(gcf,'Position',[50 50 1800 1000]);     
            clf;
            fs = 16;

            subplot(3,4,[1,2]);

            
            cla;
            x=1:hist.newt_it_number(1,sol.time_ite_save);
            y=hist.residuals(1,1:hist.newt_it_number(1,sol.time_ite_save));
            graphic = semilogy(x,y,'LineWidth',2); 
            yt=obj.logarithmic_ticks_generator(min(y),max(y));
            yticks(yt);
            ylabel("residual coupled",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[3,4]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.time_ite_save);
            y=hist.residuals(2,1:hist.newt_it_number(1,sol.time_ite_save));
            graphic = semilogy(x,y,'LineWidth',2);
            yticks(obj.logarithmic_ticks_generator(min(y),max(y)));
            ylabel("residual ps",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);

            subplot(3,4,[5,6]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.time_ite_save);
            y=hist.residuals(3,1:hist.newt_it_number(1,sol.time_ite_save));
            graphic = semilogy(x,y,'LineWidth',2);
            yticks(obj.logarithmic_ticks_generator(min(y),max(y)));
            ylabel("residual pe",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[7,8]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.time_ite_save);
            y=hist.residuals(4,1:hist.newt_it_number(1,sol.time_ite_save));
            graphic = semilogy(x,y,'LineWidth',2);
            yticks(obj.logarithmic_ticks_generator(min(y),max(y)));
            ylabel("residual ce",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[9,10]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.time_ite_save);
            y=hist.residuals(5,1:hist.newt_it_number(1,sol.time_ite_save));
            graphic = semilogy(x,y,'LineWidth',2);
            yticks(obj.logarithmic_ticks_generator(min(y),max(y)));
            ylabel("residual csn",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);
            subplot(3,4,[11,12]);
            
            cla;
            x=1:hist.newt_it_number(1,sol.time_ite_save);
            y=hist.residuals(6,1:hist.newt_it_number(1,sol.time_ite_save));
            graphic = semilogy(x,y,'LineWidth',2);
            yticks(obj.logarithmic_ticks_generator(min(y),max(y)));
            ylabel("residual csp",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Newton iteration",'FontSize',fs);

            saveas(gcf,deb.graph_folder_name+file_name);
        end

        function outcome_loc= plot_resuduals_DFN(obj,file_name)
            global hist
            global sol
            global deb
            figure(234234235);
            set(gcf,'Position',[50 50 1800 1000]);     
            clf;
            fs = 16;

            subplot(3,4,[1,2]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_time(1,1:sol.time_ite);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual coupled",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(3,4,[3,4]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_time(2,1:sol.time_ite);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual ps",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(3,4,[5,6]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_time(3,1:sol.time_ite);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual pe",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(3,4,[7,8]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_time(4,1:sol.time_ite);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual ce",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(3,4,[9,10]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_time(5,1:sol.time_ite);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual csn",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);

            subplot(3,4,[11,12]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_time(6,1:sol.time_ite);
            semilogy(x,y,'LineWidth',2);
            ylabel("residual csp",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("DFN iteration",'FontSize',fs);


            saveas(gcf,deb.graph_folder_name+file_name);
        end

        function outcome_loc= plot_resuduals_diff(obj,file_name)
            global hist
            global sol
            global deb
            figure('visible','off');
            set(gcf,'Position',[50 50 1800 1000]);     
            clf;
            fs = 16;

            subplot(4,4,[1,2]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_diff(1,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("residual coupled",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[3,4]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_diff(2,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("residual ps",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[5,6]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_diff(3,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("residual pe",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[7,8]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_diff(4,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("residual ce",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[9,10]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_diff(5,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("residual csn",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);

            subplot(4,4,[11,12]);
            cla;
            x=1:sol.time_ite;
            y=hist.residuals_diff(6,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("residual csp",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("Diff DFN iteration",'FontSize',fs);


            subplot(4,4,[13,14]);
            cla;
            x=1:sol.time_ite;
            y=hist.newt_it_number(1,1:sol.time_ite);
            plot(x,y,'LineWidth',2);
            ylabel("Newton solver iterations",'FontSize', fs);
            set(gca,'FontSize', fs);
            %legend({'$I(t)$'},'interpreter','latex');
            xlabel("iteration",'FontSize',fs);


            saveas(gcf,deb.graph_folder_name+file_name);
        end


        function  outcome_loc=plot_solid_concentration_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array)
            global fv
            global hist
            global sol
            global deb
            flag = 1;
            outcome_loc=flag;
            fs = 16;

            
            temp=size(hist.csn);
            Nn=sol.time_ite;
            temp=size(hist.csp);
            Nn2=sol.time_ite;

            
            figure('visible','off');
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
        
            


            saveas(gcf,deb.graph_folder_name+file_name);



        end


        function  outcome_loc=plot_solid_concentration_singlePart_data(obj,file_name,x,fig_nb,nnb,snb,pnb,time_array)
            global fv
            global hist
            global sol
            global deb
            flag = 1;
            outcome_loc=flag;
            fs = 16;


            temp=size(hist.csn);
            Nn=sol.time_ite;
            temp=size(hist.csp);
            Nn2=sol.time_ite;

            
            figure('visible','off');
            set(gcf,'Position',[50 50 1800 1000]);     
            
            subplot(2,2,1)
            x=sol.part_coord_n;
            y=hist.csn(1:sol.part_nb_cell+1,Nn);

            cla;
            plot(x,y,'r-+','LineWidth',2);

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
            plot(x,y,'r-+','LineWidth',2);

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
            plot(x,y,'r-+','LineWidth',2);

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
            plot(x,y,'r-+','LineWidth',2);

            if min(y) == max(y)
                ylim([min(y) min(y)+abs(min(y))*0.00001+0.000000000001]);
            else
                ylim([min(y) max(y)]);
            end

            ylabel('csp2','FontSize',fs);
            xlabel('x','FontSize',fs);
        
            title('Li concentration in pos. electrode','fontsize',fs);
            set(gca,'FontSize', fs)
        
            


            saveas(gcf,deb.graph_folder_name+file_name);


        end



    end
end
