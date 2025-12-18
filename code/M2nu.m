
function nu = M2nu(m,ecc,incl)
% --- 1. M -> E (평균 근점 이각 -> 이심 근점 이각) ---
% M = E - e*sin(E)
% 이 방정식은 E에 대해 풀 수 없으므로, 뉴턴-랩슨법(Newton-Raphson method)을 사용하여
% 수치적으로 해를 구해야 합니다.
% (입력으로 m, ecc, typeorbit이 필요하다고 가정합니다)
small =1e-10;

if ecc < small % 원궤도 계열
    if (incl < small) || (abs(incl - pi) < small)
        typeorbit = 'ce'; % Circular Equatorial
    else
        typeorbit = 'ci'; % Circular Inclined
    end
else % 비원궤도 계열
    if (incl < small) || (abs(incl - pi) < small)
        typeorbit = 'ee'; % Elliptical Equatorial
    else
        typeorbit = 'ei'; % Elliptical Inclined
    end
end


if typeorbit(1) == 'e' % ee, ei (타원 궤도)
    
    % --- 1단계: M -> E (뉴턴-랩슨법으로 E 풀기) ---
    
    % 수치해석을 위한 초기값 및 설정
    E_guess = m; % 초기 추정값 (M으로 시작하는 것이 일반적)
    max_iter = 100; % 최대 반복 횟수
    tolerance = 1e-12; % 수렴 허용 오차
    
    E = E_guess;
    for i = 1:max_iter
        f = E - ecc * sin(E) - m;  % f(E) = E - e*sin(E) - M
        f_prime = 1 - ecc * cos(E); % f'(E) = 1 - e*cos(E)
        
        delta_E = f / f_prime; % E_new = E_old - f(E)/f'(E)
        E = E - delta_E;
        
        % 수렴 조건 확인 (변화량이 허용 오차보다 작으면 종료)
        if abs(delta_E) < tolerance
            break;
        end
    end
    
    % --- 2단계: E -> nu (이심 근점 이각 -> 진 근점 이각) ---
    % 기존 코드의 역함수: nu = 2 * atan(sqrt((1+e)/(1-e)) * tan(E/2))
    nu = 2 * atan(sqrt((1+ecc)/(1-ecc)) * tan(E/2));
    
    % nu가 0 ~ 2*pi 범위에 있도록 조정
    nu = mod(nu, 2*pi);
    
    % 이 궤도에서는 arglat, truelon을 사용하지 않음
    arglat = NaN;
    truelon = NaN;

else % ci, ce (원궤도)
    % 원궤도에서는 M과 E가 정의되지 않습니다.
    % 입력된 'm' 값을 궤도상의 위치(각도)로 직접 사용합니다.
    if strcmp(typeorbit, 'ci')
        arglat = m; % m을 Argument of Latitude로 해석
        nu = NaN; 
        E = NaN;
    elseif strcmp(typeorbit, 'ce')
        truelon = m; % m을 True Longitude로 해석
        nu = NaN; 
        E = NaN;
    end
end
end