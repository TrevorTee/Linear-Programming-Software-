function Optimization = coesolver( inequals, column_num )
%   solve the inequalities and create the tableau
    coefficients = {};%每个
    optimal_row = 0; %记录objective function 的行number
    beta_col = 0;
    for index = 1:column_num %number of inequalities
        
        single_ineq = inequals{index, 1};%get the single ineq
        len = length(single_ineq);
        %初始化中间向量
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
                %i = i-1;%回退到判断下一个变量系数之前，用于for循环的增加
                
                coe_s = [coe_s;cellstr(var),cellstr(coe)] %存放系数
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
                 %存放等号右边的beta
                 while ~isstrprop(single_ineq(i),'digit')
                     i = i + 1;
                 end
                 beta = [beta,single_ineq(i:end)];
                 coe_s = [coe_s;cellstr('beta'),cellstr(beta)];
                 i = length(single_ineq) %识别完一个不等式
            elseif single_ineq(i)== '='%indentify the symbol 解决加在最后的optimal
                if single_ineq(i-1) == '<'
                    break
                else
                    coe_s = [coe_s;cellstr('beta_optimal'),cellstr('0')];
                    break;%是objective function，加一个flag
                end
            end  %if 判断的结束
        end
        if ~ismember('beta_optimal',coe_s(:,1))
            coe_s = [coe_s;cellstr(num2str(index)),cellstr('1')];
        end
        %一个不等式识别结束，把系数放到对应矩阵
        coefficients{index} = coe_s;%将系数完全放到元组当中
        
    end
    
    %得到数组的表头
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
        v_num = SIZE(1);%每一个constraints的系数个数
        for Col = 1:v_size
            for i = 1:v_num
                if strcmp(variables{Col},co{i,1})
                    Table{index_c,Col} = co{i,2}
                    break;
                elseif strcmp('beta_optimal',co{i,1})
                    if strcmp(variables{Col},'beta')
                        beta_col = Col  %识别到beta所在的列
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
    
    %系数相关矩阵Table,str到num的转换
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
        %找到indicators
        for i = 1:column_num
            if i == optimal_row
                continue;
            else
                indicator = tables(i,beta_col)/tables(i,pivot_col);
                tables(i,v_size+1) = indicator;
            end
        end
        
        %比较indicator的大小      
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
        %行变换
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

