function Optimization = coesolver( inequals, column_num )
%   solve the inequalities and create the tableau
    coefficients = {};%ÿ��
    optimal_row = 0; %��¼objective function ����number
    beta_col = 0;
    for index = 1:column_num %number of inequalities
        
        single_ineq = inequals{index, 1};%get the single ineq
        len = length(single_ineq);
        %��ʼ���м�����
        coe = [];
        var = [];
        beta = [];
        coe_s = {};
        for i = 1:len
            if isstrprop(single_ineq(i),'digit') %identify the coefficients
                coe = [coe, single_ineq(i)];
            elseif isletter(single_ineq(i))%indentify the variables
                if isempty(coe) %set 1 or -1
                    coe = '1';
                elseif (coe == '-') & (length(coe) == 1)
                    coe = '-1';
                end
                while isletter(single_ineq(i)) | isstrprop(single_ineq(i),'digit')
                    var = [var, single_ineq(i)];
                    i = i+1;
                end
                %i = i-1;%���˵��ж���һ������ϵ��֮ǰ������forѭ��������
                
                coe_s = [coe_s;cellstr(var),cellstr(coe)] %���ϵ��
                var = [];
                coe = [];
                continue;
            elseif isspace(single_ineq(i))%indentify the space
                coe = [];
                continue
            elseif single_ineq(i)== '+'%indentify the symbol
                coe = [];
                continue
            elseif single_ineq(i)== '-'%indentify the symbol
                coe = single_ineq(i);
            elseif single_ineq(i)== '<'%indentify the symbol
                 %��ŵȺ��ұߵ�beta
                 while ~isstrprop(single_ineq(i),'digit')
                     i = i + 1;
                 end
                 beta = [beta,single_ineq(i:end)];
                 coe_s = [coe_s;cellstr('beta'),cellstr(beta)];
                 i = length(single_ineq) %ʶ����һ������ʽ
            elseif single_ineq(i)== '='%indentify the symbol �����������optimal
                if single_ineq(i-1) == '<'
                    break
                else
                    coe_s = [coe_s;cellstr('beta_optimal'),cellstr('0')];
                    break;%��objective function����һ��flag
                end
            end  %if �жϵĽ���
        end
        if ~ismember('beta_optimal',coe_s(:,1))
            coe_s = [coe_s;cellstr(num2str(index)),cellstr('1')];
        end
        %һ������ʽʶ���������ϵ���ŵ���Ӧ����
        coefficients{index} = coe_s;%��ϵ����ȫ�ŵ�Ԫ�鵱��
        
    end
    
    %�õ�����ı�ͷ
    variables = {};
    for index_t = 1:column_num
        co = coefficients{index_t};
        SIZE = size(co);
        v_num = SIZE(1);
        for ind = 1:v_num
            if ismember(co{ind,1},variables)
                continue;
            else
                variables = [variables, co(ind,1)]
            end
        end
    end
    v_size = size(variables);
    v_size = v_size (2);
    
    for opt = 1:v_size
        if strcmp(variables{opt},'beta_optimal')
            if opt == v_size
                variables = [variables(1:opt-1)];
                break;
            else
                variables = [variables(1:opt-1),variables(opt+1:end)]
                break;
            end
        end
    end
    v_size = size(variables);
    v_size = v_size (2)
    
    Table = cell(column_num, v_size);
    for index_c = 1:column_num
        co = coefficients{index_c}
        SIZE = size(co);
        v_num = SIZE(1);%ÿһ��constraints��ϵ������
        for Col = 1:v_size
            for i = 1:v_num
                if strcmp(variables{Col},co{i,1})
                    Table{index_c,Col} = co{i,2}
                    break;
                elseif strcmp('beta_optimal',co{i,1})
                    if strcmp(variables{Col},'beta')
                        beta_col = Col  %ʶ��beta���ڵ���
                        Table{index_c,Col} = co{i,2};
                        optimal_row = index_c;
                    else
                        continue;
                    end
                else
                    Table{index_c,Col} = '0';
                end           
            end
        end
    end
    
    %ϵ����ؾ���Table,str��num��ת��
    table = cell(column_num, v_size);
    for in = 1:column_num
        for k = 1:v_size
            table{in,k} = str2double(Table{in,k});
        end
    end
    tables = cell2mat(table);
    %ob_f = cell2mat(table(optimal_row,:));
    while min(tables(optimal_row,:))<0
        [~, pivot_col] = min(tables(optimal_row,:));
        %�ҵ�indicators
        for i = 1:column_num
            if i == optimal_row
                continue;
            else
                indicator = tables(i,beta_col)/tables(i,pivot_col);
                tables(i,v_size+1) = indicator;
            end
        end
        
        %�Ƚ�indicator�Ĵ�С      
        [indicator,pivot_row] = max(tables(:,v_size+1));
        for t = 1:column_num
            if t == optimal_row
                continue;
            else
                if tables(t,v_size+1)<indicator
                    pivot_row = t;
                    indicator = tables(t,v_size+1)
                else continue;
                end
            end 
        end
        %�б任
        pivot = tables(pivot_row,pivot_col);
        tables(pivot_row,:) = tables(pivot_row,:)/pivot;
        for ids = 1:column_num
            if ids == pivot_row 
                continue;
            end
            
            multiplier = tables(ids,pivot_col)/pivot;
            tables(ids,:) = tables(ids,:)-tables(pivot_row,:)*multiplier;
                
        end
    end
    basic_var = {};
    basic_num = [];
    for indi = 1:v_size
        flag = 0;

        for i = 1:column_num
            if tables(i,indi) == 0
                continue;
            else
                flag = flag + 1;
            end
        end
        if flag == 1
            [~,basic_row] = max(tables(:,indi));
            if isstrprop(variables{indi},'digit')
                continue;
            end
            basic_var = [basic_var, variables(indi)];
            basic_num = [basic_num, tables(basic_row,beta_col)];
            if strcmp(variables(indi),'z')
                Max = tables(optimal_row,beta_col)/tables(optimal_row,indi);
            end
        end   
    end
    Optimization = {};
    for i = 1:length(basic_var)
       Optimization = [Optimization,cellstr([basic_var{i},'=',num2str(basic_num(i))])] 
    end
    
end

